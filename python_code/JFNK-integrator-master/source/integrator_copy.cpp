#include "points.h"
#include "integrator.h"
#include "preprocess.h"
#include <iostream>
#include <Eigen/Dense>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/tools/roots.hpp>
#include <vector>
#include <functional>
#include <set>

using namespace std;
using namespace boost::math::quadrature;
using namespace Eigen;

#define GAUSS_LEGENDRE 1
#define GAUSS_LOBATTO 2
#define GAUSS_RADAU 3
#define GAUSS_RADAU_2A 4

// 最大的JFNK迭代次数
#define N_ITER_MAX_NEWTON 50
// 牛顿-克里洛夫求解器在残差中的相对容差
#define NK_FTOL 1e-1
// 向后欧拉方法的相对容差
#define BE_TOL 1e-12
// 谱延迟校正（SDC）方法中校正强度的相对收敛准则
#define SDC_TOL 1e-14
// 最大的SDC迭代次数
#define N_ITER_MAX_SDC 500
// 最大的JFNK迭代次数
#define N_ITER_MAX_NEWTON 50
// 用于调试，自适应步长的最大步数
#define N_STEPS_MAX_ADAPTIVE 1e9
// 时间步长能增长的最大倍数
#define ADAPTIVE_SCALER_MAX 4
// 时间步长能增长的最小倍数
#define ADAPTIVE_SCALER_MIN 1.5
//
struct Params
{
    int n_nodes;            // 每个时间步的节点数
    int m;                  // 方程维度
    int node_type;          // 节点类型
    vector<double> t;       // 时间节点
    Eigen::MatrixXd S;      // 谱矩阵
    Eigen::MatrixXd S_p;    // 谱矩阵的导数
    double spectral_radius; // 谱半径

    Param(int n_nodes, int m, int node_type = GAUSS_LOBATTO)
    {
        this->n_nodes = n_nodes;
        this->m = m;
        this->node_type = node_type;
        this->t = get_nodes(n_nodes, node_type);
        this->S = spectral_matrix(t);
        this->S_p = spectral_matrix_derivative(t);
        this->spectral_radius = spectral_radius(node_type, S, S_p);
    }
}

struct jnfk_uniform_res
{
    vector<Eigen::VectorXd> y_all;          // 每个节点的解值
    vector<vector<Eigen::VectorXd>> Y_all;  // 每个时间步的迭代历史
    :vector<vector<Eigen::VectorXd>> D_all; // 每个时间步的延迟修正历史
    bool is_stiff;
    vector<double> ratios;
};

struct jnfk_step_res
{
    vector<Eigen::VectorXd> y;
    Eigen::VectorXd t;
    Eigen::VectorXd y;
    Eigen::VectorXd d;
};

// 多维牛顿迭代求解器
Eigen::VectorXd solve_root(
    const std::function<Eigen::VectorXd(Eigen::VectorXd)> &func,
    Eigen::VectorXd x0,
    double tol)
{

    int max_iter = 100;
    int iter = 0;
    Eigen::VectorXd x = x0;
    Eigen::VectorXd f = func(x);
    double res_norm = f.norm();

    while (res_norm > tol && iter < max_iter)
    {
        // 数值导数近似雅可比矩阵
        double eps = 1e-8;
        int m = x.size();
        Eigen::MatrixXd J(m, m);

        for (int j = 0; j < m; j++)
        {
            Eigen::VectorXd x_pert = x;
            x_pert(j) += eps;
            Eigen::VectorXd f_pert = func(x_pert);
            J.col(j) = (f_pert - f) / eps;
        }

        // 求解线性系统 J*dx = -f
        Eigen::VectorXd dx = J.colPivHouseholderQr().solve(-f);
        x += dx;
        f = func(x);
        res_norm = f.norm();
        iter++;
    }

    return x;
}

Eigen::VectorXd backward_euler_node(
    const function<Eigen::VectorXd(double, Eigen::VectorXd)> &f_eval,
    double t,
    double h,
    Eigen::VectorXd rhs,
    Eigen::VectorXd x0,
    double be_tol = BE_TOL)
{

    Eigen::VectorXd y = Eigen::VectorXd::Zero(rhs.size());

    if (h == 0)
    {
        y = rhs;
    }
    else
    {
        auto A = [&](Eigen::VectorXd x)
        {
            return x - h * f_eval(t, x) - rhs;
        };
        y = solve_root(A, x0, be_tol);
    }

    return y;
}

// 向后欧拉方法求解器
vector<Eigen::VectorXd> backward_euler(
    const function<Eigen::VectorXd(double, Eigen::VectorXd)> &f_eval,
    const vector<double> &t,
    const Eigen::VectorXd &y0,
    const Eigen::MatrixXd &S_p,
    double be_tol = BE_TOL)
{

    int n_nodes = t.size();
    int m = y0.size();

    // 初始化解矩阵 [n_nodes x m]
    vector<Eigen::VectorXd> y(n_nodes, y0);

    // 向后欧拉法
    for (int i = 0; i < n_nodes; i++)
    {
        double h = S_p(i, i);

        if (i == 0 && h == 0)
        {
            y[0] = y0;
        }
        else
        {

            Eigen::VectorXd rhs(m);
            if (i == 0 && h != 0)
            {
                rhs = y0;
            }
            else
            {
                rhs = y[i - 1];
            }

            y[i] = backward_euler_node(f_eval, t[i], h, rhs, y0);
        }
    }

    return y;
}

jnfk_uniform_res jfnk_uniform(
    const function<Eigrn::VectorXd(double, Eigrn::VectorXd)> f_eval,
    double t_init,
    double t_final,
    int n_steps,
    Params p,
    Eigen::VectorXd y0,
    int n_iter_max_newton = N_ITER_MAX_NEWTON,
    double be_tol = BE_TOL,
    double sdc_tol = SDC_TOL,
    bool do_print = false)
{
    assert(t_init < t_final);
    double dt = (t_final - t_init) / double(n_steps);
    Eigen::VectorXd v0 = y0;

    vector<double> t_step;
    if (p.t[p.size() - 1] != 1.0)
    {
        t_step.resize(p.size() + 1);
        std::copy(p.t.begin(), p.t.end(), t_step.begin());
        t_step.back() = 1.0;
    }
    else
    {
        t_step = p.t;
    }

    // 归一化的节点，范围在[0, 1]
    // Eigen::VectorXd t_step;
    // if (p.t(p.t.size() - 1) != 1.0) {
    //     t_step.resize(p.t.size() + 1);
    //     t_step.head(p.t.size()) = p.t;
    //     t_step(p.t.size()) = 1.0;
    // } else {
    //     t_step = p.t;
    // }

    vector<double> nodes; // 实际时间节点
    nodes.reserve(t_step.size());
    for (int i = 0; i < t_step.size(); i++)
        nodes.push_back(t_init + dt * t_step[i]);

    vector<Eigen::VectorXd> y_all;
    vector<Eigen::VectorXd> t_all;
    vector<vector<Eigen::VectorXd>> Y_all;
    vector<vector<Eigen::VectorXd>> D_all;

    // Eigen::VectorXd nodes = t_init + dt * t_step;//实际时间节点
    // Eigen::MitrixXd y_all;
    // Eigen::MatrixXd y_submatrix;
    // Eigen::VertorXd t_all;
    // vector<Eigen::MitrixXd> Y_all;
    // vector<Eigen::MitrixXd> D_all;

    y_all.push_back(y0);
    t_all.push_back(t_init);

    // y_all.row(0) = y0.transpose();
    // t_all(0) = nodes(0);

    for (int i = 0; i < n_steps; i++)
    {
        // y一次迭代Δt内的所有迭代解
        auto [y, Y, D, is_converged, is_stiff, ratios] =
            jfnk_step(f_eval, p, dt, nodes, v0, N_ITER_MAX_NEWTON, BE_TOL, SDC_TOL);

        // 验证收敛性
        assert(is_converged && "JFNK在时间步" + std::to_string(i) + "未收敛");

        for (int m = 1; m < y.rows(); m++)
        {
            y_all.push_back(y[m]);
            t_all.push_back(nodes(m));
        }
        // y_submatrix = y.bottomRows(y.rows() - 1);
        // y_all.conservativeResize(y_all.rows() + y_submatrix.rows(), Eigen::NoChange);
        // y_all.bottomRows(y_submatrix.rows()) = y_submatrix;

        Y_all.push_back(Y);
        D_all.push_back(D);
        v0 = y_all[y_all.size() - 1];

        for (int m = 0; m < t_step.size(); m++)
            nodes[m] += dt;
    }
    return {
        t_all,
        y_all,
        Y_all,
        D_all,
    };
}

tuple<vector<Eigen::VectorXd>, vector<vector<Eigen::VectorXd>>, vector<vector<Eigen::VectorXd>>, bool, bool, vector<double>>
jfnk_step(
    const function<Eigen::VectorXd(double, Eigen::VectorXd)> f_eval,
    Params p,
    double dt,
    vector<double> t,
    Eigen::VectorXd y0,
    int n_iter_max_newton = N_ITER_MAX_NEWTON,
    double be_tol = BE_TOL,
    double sdc_tol = SDC_TOL)
{
    Eigen::MatrixXd S = dt * p.S;     // 谱矩阵
    Eigen::MatrixXd S_p = dt * p.S_p; // 谱矩阵的导数
    double spectral_radius = p.spectral_radius;

    vector<Eigen::VectorXd> y_be = backward_euler(f_eval, t, y0, S_p, be_tol = be_tol);
    auto [y, Y, D, is_converged, is_stiff, ratios] =
        jfnk(f_eval, t, y0, y_be, S, S_p, spectral_radius, n_iter_max_newton, be_tol = be_tol, sdc_tol = sdc_tol);

    return {y, Y, D, is_converged, is_stiff, ratios};
}

// 计算vector<vector<Eigen::VectorXd>>结构的弗罗贝尼乌斯范数
double computeNorm(vector<vector<Eigen::VectorXd>> D)
{
    double norm_squared = 0.0;

    for (const vector<Eigen::VectorXd> i : D)
    {

        for (Eigen::VectorXd m : i)
        {
            norm_squared += m.squaredNorm();
        }
    }
    return sqrt(norm_squared);
}

tuple<vector<Eigen::VectorXd>, vector<Eigen::VectorXd>>
sdc_sweep(
    function<Eigen::VectorXd(double, const Eigen::VectorXd)> f_eval,
    vector<double> t,
    Eigen::VectorXd y0,
    vector<Eigen::VectorXd> y_old,
    Eigen::MatrixXd F,
    Eigen::MatrixXd S,
    Eigen::MatrixXd S_p,
    double be_tol = BE_TOL)
{

    int n_nodes = t.size();
    int m = y0.size();

    vector<VectorXd> y(n_nodes, VectorXd::Zero(m));
    vector<VectorXd> d(n_nodes, VectorXd::Zero(m));

    double h;
    VectorXd w;
    VectorXd rhs;
    bool do_skip;

    for (int i = 1; i < n_nodes; i++)
    {
        h = S_p(i, i);
        w = VectorXd::Zero(n_nodes);
        w(i) = h;

        do_skip = true;
        for (int m = 0; m < S_p.cols(); m++)
        {
            if (abs(S_p(i, m)) > 1e-12)
            {
                do_skip = false;
                break;
            }
        }

        if (i == 0 && do_skip)
        {
            y[0] = y0;
        }
        else
        {
            if (i == 0 && !do_skip)
            {
                rhs = y0 + ((S.row(i) - w) * F).transpose();
            }
            else
            {
                rhs = y[i - 1] + ((S.row(i) - S.row(i - 1) - w) * F).transpose();
            }

            auto [y_temp, d_temp] = sdc_node(f_eval, t[i], h, rhs, y_old[i], be_tol);
            y[i] = y_temp;
            d[i] = d_temp;
        }
    }

    return {y, d};
}

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool> sdc(
    function<Eigen::VectorXd(double, const Eigen::VectorXd)> f_eval,
    vector<double> t,
    Eigen::VectorXd y0,
    vector<Eigen::VectorXd> y_old,
    MatrixXd S,
    MatrixXd S_p,
    int n_iter_max_sdc = N_ITER_MAX_SDC,
    double be_tol = BE_TOL,
    double sdc_tol = SDC_TOL)
{

    // 获取时间节点数量和系统规模
    int n_nodes = t.size();
    int m = y0.size();

    // 存储解的历史记录
    vector<vector<Eigen::VectorXd>> Y;
    // 存储校正值的历史记录
    ector<vector<Eigen::VectorXd>> D;

    // 初始猜测解
    vector<Eigen::VectorXd> y_sdc = y_old;

    bool is_converged = false;

    // 迭代直到收敛或达到最大迭代次数
    while (!is_converged && k < n_iter_max_sdc)
    {
        // 存储当前解
        Y.push_back(y_sdc);

        Eigen::MatrixXd F(n_nodes, m);
        for (int i = 0; i < n_nodes; ++i)
        {
            // F[i] = Eigen::VectorXd::Zero(m);
            F.row(i) = f_eval(t[i], y_sdc[i]);
        }

        // 进行SDC扫描
        auto [y_sdc_new, d] = sdc_sweep(f_eval, t, y0, y_sdc, F, S, S_p, be_tol);
        is_converged = convergence_criteria(d, Y[-1], tol = sdc_tol);

        y_sdc = y_sdc_new;
        D.push_back(d);
    }

    // 返回结果：解、解的历史记录、校正值的历史记录、收敛标志
    return {y_sdc, Y, D, is_converged};
}

// JFNK初始迭代函数（用于检测系统刚性）
tuple<vector<Eigen::VectorXd>, vector<vector<Eigen::VectorXd>>, vector<vector<Eigen::VectorXd>>, bool, bool, vector<double>>
jfnk_initial(
    function<Eigen::VectorXd(double, const Eigen::VectorXd)> f_eval,
    vector<double> t,
    Eigen::VectorXd y0,
    vector<Eigen::VectorXd> y_approx,
    Eigen::MatrixXd S,
    Eigen::MatrixXd S_p,
    double spectral_radius,
    double be_tol = BE_TOL,
    double sdc_tol = SDC_TOL,
    int n_iter_max = N_ITER_MAX_SDC)
{

    bool is_stiff = false;

    // // 计数器
    // int i = 0;

    vector<Eigen::VectorXd> ratios;

    auto [y_sdc, Y_sdc, D_sdc, is_converged] =
        sdc(f_eval, t, y0, y_approx, S, S_p, 2, be_tol, sdc_tol, false);

    vector<vector<Eigen::VectorXd>> Y;
    // for (const auto &y : Y_sdc)
    // {
    //     Y.push_back(y);
    // }
    Y = Y_sdc;

    vector<Eigen::MatrixXd> D;
    // for (const auto &d : D_sdc)
    // {
    //     D.push_back(d);
    // }
    D = D_sdc;

    // JFNK初始解的迭代过程
    while (!is_converged && !is_stiff && i < n_iter_max)
    {
        vector<Eigen::VectorXd> last_D = D.back();
        vector<Eigen::VectorXd> second_last_D = D[D.size() - 2];

        // 计算它们的弗罗贝尼乌斯范数
        double norm_last = computeNorm(last_D);
        double norm_second_last = computeNorm(second_last_D);

        ratios.push_back(ratio);

        // 判断是否为刚性系统
        is_stiff = (ratio / spectral_radius) > 0.1;

        if (is_stiff)
        {
            cout << "Stiff (order convergence detected). Use JFNK." << endl;
        }
        else
        {
            cout << "Non-stiff. Use SDC" << endl;

            auto [y_sdc_new, Y_sdc_new, D_sdc_new, is_converged_new] =
                sdc(f_eval, t, y0, y_sdc, S, S_p, 1, be_tol, sdc_tol, false);

            y_sdc = y_sdc_new;
            is_converged = is_converged_new;

            // 存储近似解
            Y.push_back(Y_sdc_new[0]);

            // 存储校正值
            D.push_back(D_sdc_new[0]);
        }

        // 更新计数器
        // i++;
        return {y_sdc, Y, D, is_converged, is_stiff, ratios};
    }

    // 转换ratios为矩阵
    Eigen::MatrixXd ratios_matrix;
    if (!ratios.empty())
    {
        ratios_matrix.resize(ratios.size(), 1);
        for (size_t j = 0; j < ratios.size(); j++)
        {
            ratios_matrix(j, 0) = ratios[j](0);
        }
    }

    return {y_sdc, Y, D, is_converged, is_stiff, ratios_matrix};
}

pair<VectorXd, VectorXd> sdc_node(
    function<VectorXd(double, const VectorXd)> f_eval,
    vecotr<double> t,
    double h,
    VectorXd rhs,
    VectorXd y_old,
    double be_tol = BE_TOL)
{

    VectorXd y_new = backward_euler_node(f_eval, t, h, rhs, y_old, be_tol);
    VectorXd d = y_new-y_old;

    return {y_new, d};
}

// JFNK迭代函数（无雅可比牛顿迭代）
tuple<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>, bool>
jfnk_iterations(
    function<Eigen::VectorXd(double, Eigen::VectorXd)> f_eval,
    vector<double> t,
    Eigen::VectorXd y0,
    vector<Eigen::VectorXd> y_init,
    Eigen::MatrixXd S,
    Eigen::MatrixXd S_p,
    int n_iter_max_newton = N_ITER_MAX_NEWTON,
    double be_tol = BE_TOL,
    double sdc_tol = SDC_TOL)
{

    vector<Eigen::VectorXd> y = y_init;

    int k_sdc = t.size() + 1;

    int m = y0.size();

    vector<vector<Eigen::VectorXd>> D_list, Y_list;
    vector<vector<Eigen::VectorXd>>

    int k = 0;
    bool is_converged = false;

    // 牛顿迭代过程
    while (!is_converged && k < n_iter_max_newton)
    {
        // 进行SDC迭代
        auto [y_new, Y, D, is_converged_new] =
            sdc(f_eval, t, y0, y, S, S_p, k_sdc, be_tol, sdc_tol);

        y = y_new;

        if (is_converged_new)
        {
            Y_list.insert(Y_list.end(), Y.begin(), Y.end());
            D_list.insert(D_list.end(), D.begin(), D.end());
            is_converged = true;
        }
        else
        {
            Y_list.insert(Y_list.end(), Y.begin(), prev(Y.end()));
            D_list.insert(D_list.end(), D.begin(), prev(D.end()));

            y = Y.back();

        }

        for (int i = 0; i < m; i++) {

        }

        // 更新迭代次数
        k++;
    }

    // 展平历史记录
    std::vector<Eigen::MatrixXd> Y_flat, D_flat;
    for (const auto &Y_sub : Y_list)
    {
        Y_flat.insert(Y_flat.end(), Y_sub.begin(), Y_sub.end());
    }
    for (const auto &D_sub : D_list)
    {
        D_flat.insert(D_flat.end(), D_sub.begin(), D_sub.end());
    }

    return {y, Y_flat, D_flat, is_converged};
}

tuple<vector<Eigen::VectorXd>, vector<vector<Eigen::VectorXd>>, vector<vector<Eigen::VectorXd>>, bool, bool, vector<double>>
jfnk(
    const function<Eigen::VectorXd(double, Eigen::VectorXd)> f_eval,
    vector<double> t,
    Eigen::VectorXd y0,
    vector<Eigen::VectorXd> y_approx,
    Eigen::MatrixXd S,
    Eigen::MatrixXd S_p,
    double spectral_radius,
    int n_iter_max_newton = N_ITER_MAX_NEWTON,
    double be_tol = BE_TOL,
    double sdc_tol = SDC_TOL)
{

    auto [y_init, Y_init, D_init, is_converged, is_stiff, ratios] = jfnk_initial(f_eval, t, y0, y_approx, S, S_p, spectral_radius, n_iter_max_newton, be_tol, sdc_tol);

    if (is_stiff)
    {
        auto [y_final, Y_final, D_final, is_converged_final] =
            jfnk_iterations(f_eval, t, y0, y_init, S, S_p, n_iter_max_newton, be_tol, sdc_tol);

        return {y_final, Y_final, D_final, is_converged_final, true, ratios_initial};
    }
    else
    {
        return {y_init, Y_init, D_init, is_converged, is_stiff, ratios};
    }
}
