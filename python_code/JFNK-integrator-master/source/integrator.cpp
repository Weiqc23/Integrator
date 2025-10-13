#include "points.h"
#include "integrator.h"
#include "preprocess.h"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/tools/roots.hpp>
#include <vector>
#include <functional>
#include <set>
#include <chrono>

using namespace std;
using namespace boost::math::quadrature;
using namespace Eigen;

#define GAUSS_LEGENDRE 1
#define GAUSS_LOBATTO 2
#define GAUSS_RADAU 3
#define GAUSS_RADAU_2A 4

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
Params::Params(int n_nodes, int m, int node_type)
{
    this->n_nodes = n_nodes;
    this->m = m;
    this->node_type = node_type;

    // 用 points.h 里的函数生成节点和谱矩阵
    this->t = get_nodes(n_nodes, node_type);
    this->S = spectral_matrix(t);
    this->S_p = backward_euler_matrix(t);
    this->spectral_radius = ::spectral_radius(node_type, S, S_p);
}

// hybr
// VectorXd root(
//     function<VectorXd(const VectorXd &)> f,
//     const VectorXd &x0,
//     double tol,
//     int max_iter,
//     bool do_print)
// {
//     auto start = chrono::high_resolution_clock::now();

//     int m = x0.size();
//     VectorXd x = x0;
//     VectorXd Fx = f(x);

//     if (Fx.norm() < tol)
//         return x;

//     // ---------- 初始 Jacobian（差分法） ----------
//     MatrixXd J(m, m);
//     const double REL_DELTA = 1e-6;
//     for (int k = 0; k < m; ++k)
//     {
//         VectorXd x_partial = x;
//         double eps = REL_DELTA * max(1.0, abs(x(k)));
//         x_partial(k) += eps;
//         VectorXd F_partial = f(x_partial);
//         J.col(k) = (F_partial - Fx) / eps;
//     }

//     // ---------- 迭代 ----------
//     for (int iter = 0; iter < max_iter; ++iter)
//     {
//         // 解 J * delta = -F
//         VectorXd delta = J.colPivHouseholderQr().solve(-Fx);

//         // 如果步长太大，缩放一下（信赖域近似）
//         double step_norm = delta.norm();
//         if (step_norm > 1.0)
//         {
//             delta *= (1.0 / step_norm);
//         }

//         // 更新解
//         VectorXd x_new = x + delta;
//         VectorXd Fx_new = f(x_new);

//         // 收敛判据
//         if (Fx_new.norm() < tol || delta.norm() < tol)
//         {
//             x = x_new;
//             break;
//         }

//         // ---------- Broyden rank-1 更新 ----------
//         VectorXd s = x_new - x;
//         VectorXd y = Fx_new - Fx;
//         double denom = s.dot(s);
//         if (denom > 1e-14)
//         {
//             J += ((y - J * s) * s.transpose()) / denom;
//         }

//         // 准备下一轮
//         x = x_new;
//         Fx = Fx_new;
//     }

//     if (do_print)
//     {
//         auto end = chrono::high_resolution_clock::now();
//         chrono::duration<double> duration = end - start;
//         cout << "root(hybr) run time: " << duration.count() << " s" << endl;
//     }

//     return x;
// }
// VectorXd root(
//     function<VectorXd(const VectorXd &)> f,
//     const VectorXd &x0,
//     double tol,
//     int max_iter,
//     bool do_print)
// {
//     auto start = std::chrono::high_resolution_clock::now();

//     int m = x0.size();
//     VectorXd x = x0;
//     VectorXd Fx = f(x);

//     if (Fx.norm() < tol) return x;

//     // 初始化 Jacobian 用差分法
//     MatrixXd J(m, m);
//     const double REL_DELTA = 1e-6;
//     for (int k = 0; k < m; ++k)
//     {
//         VectorXd x_partial = x;
//         double eps = REL_DELTA * std::max(1.0, std::abs(x(k)));
//         x_partial(k) += eps;
//         VectorXd F_partial = f(x_partial);
//         J.col(k) = (F_partial - Fx) / eps;
//     }

//     for (int iter = 0; iter < max_iter; ++iter)
//     {
//         // 解 J * delta = -Fx
//         VectorXd delta = J.colPivHouseholderQr().solve(-Fx);

//         if (delta.norm() < tol) break;

//         VectorXd x_new = x + delta;
//         VectorXd Fx_new = f(x_new);

//         if (Fx_new.norm() < tol)
//         {
//             x = x_new;
//             break;
//         }

//         // Broyden 更新
//         VectorXd s = x_new - x;
//         VectorXd y = Fx_new - Fx;
//         double denom = s.dot(s);
//         if (denom > 1e-14)
//         {
//             J += ((y - J * s) * s.transpose()) / denom;
//         }

//         x = x_new;
//         Fx = Fx_new;
//     }

//     if (do_print)
//     {
//         auto end = std::chrono::high_resolution_clock::now();
//         std::chrono::duration<double> duration = end - start;
//         cout << "root run time: " << duration.count() << " s" << endl;
//     }

//     return x;
// }

// VectorXd root(
//     function<VectorXd(const VectorXd &)> f,
//     const VectorXd &x0,
//     double tol,
//     int max_iter,
//     bool do_print)
// {

//     const double EPS_MACHINE = std::numeric_limits<double>::epsilon(); // ~2.22e-16
//     // 选择更稳健的差分步长常量（你可以调为1e-7/1e-6等做测试）
//     const double REL_DELTA = 1e-6;
//     auto start = std::chrono::high_resolution_clock::now();

//     VectorXd x = x0;
//     int m = x.size();

//     // 预分配，避免循环中重复分配
//     VectorXd Fx(m);
//     VectorXd x_partial(m);
//     MatrixXd J = MatrixXd::Zero(m, m);

//     for (int iter = 0; iter < max_iter; ++iter)
//     {
//         // 计算 Fx
//         Fx = f(x);

//         // 检查收敛
//         if (Fx.norm() < tol)
//         {
//             // 计算并打印运行时间
//             if (do_print)
//             {
//                 auto end = std::chrono::high_resolution_clock::now();
//                 std::chrono::duration<double> duration = end - start;
//                 cout << "root run time: " << duration.count() << " s" << endl;
//             }

//             return x;
//         }

//         // 计算差分雅可比（前向差分），使用相对步长
//         // eps_k = REL_DELTA * max(1, |x_k|)
//         for (int k = 0; k < m; ++k)
//         {
//             x_partial = x;
//             double xi = x(k);
//             double eps = REL_DELTA * std::max(1.0, std::abs(xi));
//             x_partial(k) += eps;
//             VectorXd F_partial = f(x_partial);
//             // 列赋值
//             J.col(k) = (F_partial - Fx) / eps;
//         }

//         // solve J * delta = -Fx
//         // 使用 ColPivHouseholderQr (通常比 fullPivLU 更快)
//         VectorXd delta;
//         // 如果 J 可能奇异，可以考虑 colPivHouseholderQr().solve(...)
//         delta = J.colPivHouseholderQr().solve(-Fx);

//         // 更新 x（可做简单阻尼以防 divergence）
//         // 如果 delta 很大，做阻尼（可选）
//         double max_delta = delta.lpNorm<Infinity>();
//         if (std::isfinite(max_delta) && max_delta > 1e2)
//         {
//             delta *= (1e2 / max_delta);
//         }

//         x += delta;

//         if (delta.norm() < tol)
//         {
//             // 计算并打印运行时间
//             if (do_print)
//             {
//                 auto end = std::chrono::high_resolution_clock::now();
//                 std::chrono::duration<double> duration = end - start;
//                 cout << "root run time: " << duration.count() << " s" << endl;
//             }
//             return x;
//         }
//     }

//     // 若未收敛，返回当前值
//     // 计算并打印运行时间
//     if (do_print)
//     {
//         auto end = std::chrono::high_resolution_clock::now();
//         std::chrono::duration<double> duration = end - start;
//         cout << "root run time: " << duration.count() << " s" << endl;
//     }
//     return x;
// }

// 手写求解器f(x)=0
VectorXd root(
    function<VectorXd(const VectorXd &)> f,
    const VectorXd &x0,
    double tol,
    int max_iter,
    bool do_print)
{

    VectorXd x = x0;
    int m = x.size();

    for (int i = 0; i < max_iter; i++)
    {
        VectorXd Fx = f(x);
        // cout << "root i: " << i << endl;
        // cout << "x: " << x << endl;
        // cout << "Fx: " << Fx << endl;
        // cout << "Fx.norm(): " << double(Fx.norm()) << endl;
        if (Fx.norm() < tol)
        {
            return x;
        }

        MatrixXd J(m, m);
        double eps = 1e-8;
        for (int k = 0; k < m; k++)
        {
            VectorXd x_partial = x;
            x_partial(k) += eps;
            VectorXd F_partial = f(x_partial);
            J.col(k) = (F_partial - Fx) / eps;
        }

        VectorXd delta = J.fullPivLu().solve(-Fx);
        x += delta;

        // cout << "delta.norm(): " << double(delta.norm()) << endl;
        if (delta.norm() < tol)
        {
            return x;
        }
    }

    return x;
}

VectorXd backward_euler_node(
    function<VectorXd(double, const VectorXd &)> f_eval,
    double t,
    double h,
    const VectorXd &rhs,
    const VectorXd &x0,
    double be_tol)
{

    VectorXd y = VectorXd::Zero(rhs.size());

    if (h == 0.0)
    {
        return rhs;
    }

    auto F = [&](const VectorXd &x)
    {
        VectorXd val = x - h * f_eval(t, x) - rhs;
        // cout << "inside F: f_eval=" << val.transpose() << endl;
        return val.eval();
    };
    // cout << "t: " << t << endl;
    // cout << "h: " << h << endl;
    // cout << "rhs.size(): " << rhs.size() << endl;
    // cout << "rhs: " << rhs << endl;
    // cout << "rhs full precision: " << std::setprecision(16) << rhs << endl;

    // cout << "f_eval: " << f_eval(t, x0) << endl;
    // cout << "F的类型ID：" << typeid(F).name() << endl;
    // cout << "F(x0): " << F(x0) << endl;
    // cout << "F(x0)(0) full precision: " << std::setprecision(16) << F(x0)(0) << endl;
    y = root(F, x0, be_tol);
    // cout << "root函数返回结果: " << endl;
    // cout << "y.size(): " << y.size() << endl;
    // cout << y << endl;
    // cout << "x0.size(): " << x0.size() << endl;
    // cout << "x0: " << x0 << endl;
    return y;
}

MatrixXd backward_euler(
    function<VectorXd(double, const VectorXd &)> f_eval,
    const VectorXd &t,
    const VectorXd &y0,
    const MatrixXd &S_p,
    double be_tol,
    bool do_print)
{
    // 计算并打印运行时间
    auto start = std::chrono::high_resolution_clock::now();

    int n_nodes = t.size();
    int m = y0.size();

    MatrixXd y = MatrixXd::Zero(n_nodes, m);
    // cout << "start backward_euler" << endl;
    for (int i = 0; i < n_nodes; i++)
    {
        double h = S_p(i, i);

        // cout << endl << "i: " << i << endl;
        // cout << "h: " << h << endl;
        if (i == 0 && h == 0)
        {
            y.row(0) = y0.transpose();
        }
        else
        {
            VectorXd rhs(m);
            if (i == 0 && h != 0)
            {
                rhs = y0;
            }
            else
            {
                rhs = y.row(i - 1).transpose();
            }
            y.row(i) = (backward_euler_node(f_eval, t(i), h, rhs, rhs)).transpose();
        }

        // cout << "y.row(i): " << i << endl;
        // cout << y.row(i) << endl;
    }
    // 计算并打印运行时间
    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "backward_euler run time: " << duration.count() << " s" << endl;
    }
    return y;
}

double adjust_scaler(double x, double x_min = ADAPTIVE_SCALER_MIN, double x_max = ADAPTIVE_SCALER_MAX)
{
    if (x > x_min)
    {
        return min(x, x_max);
    }
    else if (x >= 1 && x <= x_min)
    {
        return 1.0;
    }

    return x;
}

double scaler_lobatto(double aerr, int n_nodes, int k, double tol)
{

    int p = 2 * n_nodes - 1;
    // cout << "p: " << p << endl;
    double x = tol * 1.0 * (1 - pow(double(1.0 / k), double(p))) / double(aerr);
    // cout << "x: " << x << endl;
    double scaler = pow(x, double(1.0 / p));
    // cout << "scaler_lobatto: " << scaler << endl;
    return scaler;
}

double step_size_scaler(const VectorXd &y, const VectorXd &ysteps, int n_nodes, int k, double tol, int node_type)
{

    assert(node_type == GAUSS_LOBATTO);

    double aerr = (y - ysteps).lpNorm<Infinity>();
    // cout << "aerr: " << aerr << endl;
    double scaler;

    if (node_type == GAUSS_LOBATTO)
    {
        scaler = scaler_lobatto(aerr, n_nodes, k, tol);
        // cout << "scaler: " << scaler << endl;
    }

    return scaler;
}

double update_step_size(double dt, double t0, double t_final)
{
    double dt_new = dt;

    if ((t0 + dt) > t_final)
    {
        dt_new = t_final - t0;
    }

    return dt_new;
}

double relative_norm(const MatrixXd &top, const MatrixXd &bot)
{
    double top_norm = top.norm(); // Frobenius 范数
    double bot_norm = bot.norm();

    double ratio = top_norm / bot_norm;

    if (bot_norm == 0.0)
    {
        return (top_norm == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
    }

    return top_norm / bot_norm;
}
bool convergence_criteria(const MatrixXd &d, const MatrixXd &y, double tol)
{
    double value = relative_norm(d, y);
    bool is_converged = (value <= tol);
    return is_converged;
}

bool stopping_criteria(double t_final, double t0)
{
    double T_SMALL = 1e-15;
    double EPS = 1e-20;
    double t = t0 + EPS;
    bool do_stop;

    if (t_final >= T_SMALL || t_final == 0)
    {
        do_stop = (t >= t_final);
    }
    else
    {
        do_stop = ((t_final - t) / t_final <= 1e-5);
    }

    return do_stop;
}

pair<VectorXd, VectorXd> sdc_node(
    function<VectorXd(double, const VectorXd &)> f_eval,
    double t,
    double h,
    const VectorXd &rhs,
    const VectorXd &y_old,
    double be_tol,
    bool do_print)
{

    auto start = std::chrono::high_resolution_clock::now();

    VectorXd y_new = backward_euler_node(f_eval, t, h, rhs, y_old, be_tol);
    VectorXd d = y_new - y_old;

    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "sdc node run time: " << duration.count() << " s" << endl;
    }
    return {y_new, d};
}

tuple<MatrixXd, MatrixXd>
sdc_sweep(
    function<VectorXd(double, const VectorXd &)> f_eval,
    const VectorXd &t,
    const VectorXd &y0,
    const MatrixXd &y_old,
    const MatrixXd &F,
    const MatrixXd &S,
    const MatrixXd &S_p,
    double be_tol,
    bool do_print)
{

    auto start = std::chrono::high_resolution_clock::now();

    int n_nodes = t.size();
    int m = y0.size();

    MatrixXd y = MatrixXd::Zero(n_nodes, m);
    MatrixXd d = MatrixXd::Zero(n_nodes, m);

    double h;
    VectorXd w;
    VectorXd rhs;
    bool do_skip;
    // cout << "sdc sweep run" << endl;
    // cout << "y0: " << endl;
    // cout << y0 << endl;
    // cout << "S: " << endl;
    // cout << S << endl;
    for (int i = 0; i < n_nodes; i++)
    {
        // cout << "i: " << i << endl;
        h = S_p(i, i);
        w = VectorXd::Zero(n_nodes);
        w(i) = h;

        do_skip = (S_p.row(i).isZero(1e-15));

        if (i == 0 && do_skip)
        {
            y.row(i) = y0.transpose();
        }
        else
        {
            // cout << "y0: " << endl;
            // cout << y0 << endl;
            // cout << "S: " << endl;
            // cout << S << endl;
            // cout << "do_skip: " << do_skip << endl;
            if (i == 0 && !do_skip)
            {
                rhs = y0 + ((S.row(i) - w.transpose()) * F).transpose();
            }
            else
            {
                rhs = y.row(i - 1).transpose() + ((S.row(i) - S.row(i - 1) - w.transpose()) * F).transpose();
                // rhs = y.row(i - 1).transpose() + ((S.row(i) - S.row(i - 1) - w.transpose()) * F);
            }

            // cout << "rhs: " << endl;
            // cout << rhs << endl;
            // cout << "w: " << endl;
            // cout << w << endl;
            // cout << "F: " << endl;
            // cout << F << endl;
            // cout << "val: " << endl;
            // cout << ((S.row(i) - S.row(i - 1) - w.transpose()) * F).transpose() << endl;
            auto [y_temp, d_temp] = sdc_node(f_eval, t(i), h, rhs, y_old.row(i).transpose(), be_tol);
            y.row(i) = y_temp.transpose();
            d.row(i) = d_temp.transpose();
        }

        // cout << "y: " << endl;
        // cout << y.row(i) << endl;
    }

    // cout << "sdc_sweep y: " << endl;
    // cout << y << endl;
    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "sdc sweep run time: " << duration.count() << " s" << endl;
    }
    return {y, d};
}

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool>
sdc(
    function<VectorXd(double, const VectorXd &)> f_eval,
    const VectorXd &t,
    const VectorXd &y0,
    const MatrixXd &y_old,
    const MatrixXd &S,
    const MatrixXd &S_p,
    int n_iter_max_sdc,
    double be_tol,
    double sdc_tol,
    bool do_print)
{

    auto start = std::chrono::high_resolution_clock::now();

    int n_nodes = t.size();
    int m = y0.size();
    MatrixXd F = MatrixXd::Zero(n_nodes, m);

    vector<MatrixXd> Y;
    vector<MatrixXd> D;

    MatrixXd y_sdc = y_old;
    MatrixXd d;

    bool is_converged = false;

    int k = 0;

    while (!is_converged && k < n_iter_max_sdc)
    {
        Y.push_back(y_sdc);

        F = MatrixXd::Zero(n_nodes, m);
        for (int i = 0; i < n_nodes; i++)
        {
            F.row(i) = f_eval(t[i], y_sdc.row(i).transpose()).transpose();
            // if (n_iter_max_sdc > 2)
            // {
            // cout << endl
            //      << "build F: i: " << i << endl;
            // cout << "t[i]: " << t[i] << endl;
            // cout << "y_sdc[i]: " << endl;
            // cout << y_sdc.row(i) << endl;
            // cout << "F[i]: " << endl;
            // cout << F.row(i) << endl;
            // }
        }

        // if (n_iter_max_sdc > 2)
        // {
        //     cout << "\n\nstart sweep k: " << k << endl;
        //     tie(y_sdc, d) = sdc_sweep(f_eval, t, y0, y_sdc, F, S, S_p, be_tol, true);
        //     cout << "\ny_sdc: " << endl;
        //     cout << y_sdc << endl;
        //     cout << endl;
        // }
        // else
        // {
        tie(y_sdc, d) = sdc_sweep(f_eval, t, y0, y_sdc, F, S, S_p, be_tol);
        // }
        is_converged = convergence_criteria(d, Y.back(), sdc_tol);

        D.push_back(d);
        k++;
    }
    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "sdc run time: " << duration.count() << " s" << endl;
    }
    return {y_sdc, Y, D, is_converged};
}

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool, bool, VectorXd>
jfnk_initial(
    function<VectorXd(double, const VectorXd &)> f_eval,
    const VectorXd &t,        // 时间节点
    const VectorXd &y0,       // 初始条件
    const MatrixXd &y_approx, // 初始近似解
    const MatrixXd &S,        // 谱积分矩阵
    const MatrixXd &S_p,      // 向后欧拉积分矩阵
    double spectral_radius,   // 谱半径
    double be_tol,
    double sdc_tol,
    int n_iter_max,
    bool do_print)
{

    auto start = std::chrono::high_resolution_clock::now();

    bool is_stiff = false;
    vector<double> ratios_vec;

    MatrixXd y_sdc;
    vector<MatrixXd> Y_sdc;
    vector<MatrixXd> D_sdc;
    bool is_converged = false;

    // cout << "jfnk_initial i: " << -2 << endl;
    // cout << "y_approx" << endl;
    // cout << y_approx << endl;

    tie(y_sdc, Y_sdc, D_sdc, is_converged) =
        sdc(f_eval, t, y0, y_approx, S, S_p, 2, be_tol, sdc_tol);

    // cout << "jfnk_initial i: " << -1 << endl;
    // cout << "y_sdc" << endl;
    // cout << y_sdc << endl;

    vector<MatrixXd> Y(Y_sdc.begin(), Y_sdc.end());
    vector<MatrixXd> D(D_sdc.begin(), D_sdc.end());

    // if (is_converged)
    // {
    //     // ratios空
    // }
    int i = 0;
    // cout << "sdc run" << endl;
    while ((!is_converged) && (!is_stiff) && (i < n_iter_max))
    {
        double norm_last = D.back().norm();
        double norm_prev = D[D.size() - 2].norm();
        double ratio = norm_last / norm_prev;
        ratios_vec.push_back(ratio);

        is_stiff = (ratio / spectral_radius) > 0.1;

        if (is_stiff)
        {
            // cout << "is stiff" << endl;
        }
        else
        {
            // cout << "not stiff" << endl;

            tie(y_sdc, Y_sdc, D_sdc, is_converged) =
                sdc(f_eval, t, y0, y_sdc, S, S_p, 1, be_tol, sdc_tol);

            Y.push_back(Y_sdc[0]);
            D.push_back(D_sdc[0]);
        }
        // cout << "jfnk_initial i: " << i << endl;
        // cout << "y_sdc" << endl;
        // cout << y_sdc << endl;
        i++;
    }

    VectorXd ratios(ratios_vec.size());
    for (int k = 0; k < (int)ratios_vec.size(); ++k)
    {
        ratios(k) = ratios_vec[k];
    }

    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "jfnk_initial run time: " << duration.count() << " s" << endl;
    }
    return {y_sdc, Y, D, is_converged, is_stiff, ratios};
}

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool>
jfnk_iterations(
    function<VectorXd(double, const VectorXd &)> f_eval,
    const VectorXd &t,
    const VectorXd &y0,
    const MatrixXd &y_init,
    const MatrixXd &S,
    const MatrixXd &S_p,
    int n_iter_max_newton,
    double be_tol,
    double sdc_tol,
    bool do_print)
{

    auto start = std::chrono::high_resolution_clock::now();

    MatrixXd y = y_init;
    // cout << "y_init: " << endl;
    // cout << y_init << endl;
    int k_sdc = t.size() + 1;
    int m = y0.size();

    vector<MatrixXd> Y, D;
    vector<MatrixXd> D_list, Y_list;
    bool is_converged = false;

    int k = 0;

    while (!is_converged && k < n_iter_max_newton)
    {
        tie(y, Y, D, is_converged) = sdc(f_eval, t, y0, y, S, S_p, k_sdc, be_tol, sdc_tol);
        // cout << "k: " << k << " / " << n_iter_max_newton << endl;
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        // cout << "jfnk run time: " << duration.count() << " s" << endl;
        // cout << "is_converged: " << is_converged << endl;
        // cout << "y: " << endl;
        // cout << y << endl;
        if (is_converged)
        {
            Y_list.insert(Y_list.end(), Y.begin(), Y.end());
            D_list.insert(D_list.end(), D.begin(), D.end());
        }
        else
        {
            // auto is_converged_start = std::chrono::high_resolution_clock::now();
            // cout << "\n!is_not_converged" << endl;
            Y_list.insert(Y_list.end(), Y.begin(), Y.end() - 1);
            D_list.insert(D_list.end(), D.begin(), D.end() - 1);
            y = Y.back();

            MatrixXd yT = y.transpose();
            vector<MatrixXd> DT;
            for (int i = 0; i < D.size(); i++)
            {
                DT.push_back(D[i].transpose());
            }

            vector<MatrixXd> A;
            for (int i = 0; i < k_sdc - 1; i++)
            {
                A.push_back(DT[i + 1] - DT[i]);
            }

            vector<MatrixXd> M;
            for (int i = 0; i < k_sdc - 1; i++)
            {
                M.push_back(DT[i]);
            }

            for (int i = 0; i < m; i++)
            {
                MatrixXd B(A.size(), A[0].cols());
                for (int k = 0; k < A.size(); k++)
                {
                    B.row(k) = A[k].row(i);
                }
                MatrixXd BT = B.transpose();
                // B = B.transpose();
                VectorXd rhs = -DT.back().row(i).transpose();

                // 最小二乘解，第i个分量的c
                VectorXd c = BT.colPivHouseholderQr().solve(rhs);

                MatrixXd V(M.size(), M[0].cols());
                for (int k = 0; k < M.size(); k++)
                {
                    V.row(k) = M[k].row(i);
                }

                VectorXd dy = V.transpose() * c;
                // 更新 y的第i个分量
                yT.row(i) += dy.transpose();
            }

            y = yT.transpose();
            k++;
            // auto is_converged_end = std::chrono::high_resolution_clock::now();
            // std::chrono::duration<double> is_converged_duration = is_converged_end - is_converged_start;
            // cout << "is_converged run time: " << is_converged_duration.count() << " s" << endl;
            // cout << "y: " << endl;
            // cout << y << endl;
            // cout << endl;
        }
    }
    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "jfnk_iterations run time: " << duration.count() << " s" << endl;
    }

    return {y, Y_list, D_list, is_converged};
}

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool, bool, VectorXd>
jfnk(
    function<VectorXd(double, const VectorXd &)> f_eval, // f(t,y) 函数
    const VectorXd &t,                                   // 时间节点
    const VectorXd &y0,                                  // 初始条件
    const MatrixXd &y_approx,                            // 初始近似解
    const MatrixXd &S,                                   // 谱积分矩阵
    const MatrixXd &S_p,                                 // 预条件积分矩阵
    double spectral_radius,                              // 谱半径
    int n_iter_max_newton,
    double be_tol,
    double sdc_tol,
    bool do_print)
{

    auto start = std::chrono::high_resolution_clock::now();

    MatrixXd y_init;
    vector<MatrixXd> Y_init;
    vector<MatrixXd> D_init;
    bool is_converged = false;
    bool is_stiff = false;
    VectorXd ratios;

    // 初始 SDC 迭代
    tie(y_init, Y_init, D_init, is_converged, is_stiff, ratios) =
        jfnk_initial(f_eval, t, y0, y_approx, S, S_p, spectral_radius, be_tol, sdc_tol);

    // cout << "jfnk_initial" << endl;
    // cout << "t.size(): " << t.size() << endl;
    // cout << t << endl;
    // cout << "y0.size(): " << y0.size() << endl;
    // cout << y0 << endl;
    // cout << "y_init.size(): " << y_init.size() << endl;
    // cout << y_init << endl;
    vector<MatrixXd> Y, D;
    vector<MatrixXd> Y_all, D_all;
    MatrixXd y;

    if (is_converged)
    {
        y = y_init;
    }

    if ((!is_converged) && is_stiff)
    {
        tie(y, Y, D, is_converged) =
            jfnk_iterations(f_eval, t, y0, y_init, S, S_p, n_iter_max_newton, be_tol, sdc_tol);
    }

    vector<MatrixXd> Y_new;
    vector<MatrixXd> D_new;
    Y_new.insert(Y_new.end(), Y_init.begin(), Y_init.end());
    Y_new.insert(Y_new.end(), Y.begin(), Y.end());
    Y = Y_new;

    D_new.insert(D_new.end(), D_init.begin(), D_init.end());
    D_new.insert(D_new.end(), D.begin(), D.end());
    D = D_new;

    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "jfnk run time: " << duration.count() << " s" << endl;
    }
    return {y, Y, D, is_converged, is_stiff, ratios};
}

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool, bool, VectorXd>
jfnk_step(
    function<VectorXd(double, const VectorXd &)> f_eval,
    const Params &p,
    double dt,
    const VectorXd &t,
    const VectorXd &y0,
    int n_iter_max_newton,
    double be_tol,
    double sdc_tol,
    bool do_print)
{
    auto start = std::chrono::high_resolution_clock::now();

    MatrixXd S = dt * p.S;
    MatrixXd S_p = dt * p.S_p;
    double spectral_radius = p.spectral_radius;

    MatrixXd y_be = backward_euler(f_eval, t, y0, S_p, be_tol);
    // cout << "y0: " << endl;
    // cout << y0 << endl;
    // cout << "y_be: " << endl;
    // cout << y_be << endl;
    // cout << "S_p" << endl;
    // cout << S_p << endl;
    // cout << "t.size(): " << t.size() << endl;
    // cout << t << endl;
    // cout << "y0.size(): " << y0.size() << endl;
    // cout << y0 << endl;
    // cout << "backward_euler返回结果: " << endl;
    // cout << "y.size(): " << y_be.size() << endl;
    // cout << y_be.cast<double>() << endl;

    MatrixXd y;
    vector<MatrixXd> Y;
    vector<MatrixXd> D;
    bool is_converged;
    bool is_stiff;
    VectorXd ratios;

    tie(y, Y, D, is_converged, is_stiff, ratios) =
        jfnk(f_eval, t, y0, y_be, S, S_p, spectral_radius, n_iter_max_newton, be_tol = be_tol, sdc_tol = sdc_tol);

    // cout << "t.size(): " << t.size() << endl;
    // cout << t << endl;
    // cout << "y.size(): " << y.size() << endl;
    // cout << y << endl;
    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "jfnk_step run time: " << duration.count() << " s" << endl;
    }
    return {y, Y, D, is_converged, is_stiff, ratios};
}

tuple<VectorXd, MatrixXd, vector<vector<MatrixXd>>, vector<vector<MatrixXd>>>
jfnk_uniform(
    function<VectorXd(double, const VectorXd &)> f_eval,
    double t_init,
    double t_final,
    int n_steps,
    const Params &p,
    const VectorXd &y0,
    int n_iter_max_newton,
    double be_tol,
    double sdc_tol,
    bool do_print)
{
    auto start = std::chrono::high_resolution_clock::now();
    assert(t_init < t_final);
    double dt = (t_final - t_init) / double(n_steps);

    VectorXd v0 = y0;
    VectorXd t_step;

    if (p.t(p.t.size() - 1) != 1.0)
    {
        t_step.resize(p.t.size() + 1);
        t_step.head(p.t.size()) = p.t;
        t_step(p.t.size()) = 1.0;
    }
    else
    {
        t_step = p.t;
    }

    VectorXd nodes(t_step.size());
    nodes = dt * t_step + VectorXd::Ones(t_step.size()) * t_init;

    MatrixXd y_all;
    VectorXd t_all;
    MatrixXd y_sub;
    VectorXd t_sub;
    vector<vector<MatrixXd>> Y_all;
    vector<vector<MatrixXd>> D_all;

    y_all.conservativeResize(1, y0.size());
    y_all.row(0) = y0.transpose();
    t_all.conservativeResize(1, 1);
    t_all(0) = t_init;

    for (int i = 0; i < n_steps; i++)
    {
        MatrixXd y;
        vector<MatrixXd> Y;
        vector<MatrixXd> D;
        bool is_converged;
        bool is_stiff;
        VectorXd ratios;
        // if (n_steps>1) {
        //     cout << "n_steps k: " << i << endl;
        //     cout << "v0: " << endl;
        //     cout << v0 << endl;
        //     cout << "dt: " << dt << endl;
        //     cout << "nodes: " << endl;
        //     cout << nodes << endl;
        // }
        tie(y, Y, D, is_converged, is_stiff, ratios) =
            jfnk_step(f_eval, p, dt, nodes, v0, n_iter_max_newton, be_tol, sdc_tol);

        // if (n_steps>1) {
        //     cout << "k: " << i << endl;
        //     cout << "y: " << endl;
        //     cout << y << endl;
        // }

        assert(is_converged && (std::string("JFNK在时间步") + std::to_string(i) + "未收敛").c_str());

        y_sub = y.block(1, 0, y.rows() - 1, y.cols());
        y_all.conservativeResize(y_all.rows() + y_sub.rows(), NoChange);
        y_all.bottomRows(y_sub.rows()) = y_sub;

        t_sub = nodes.block(1, 0, nodes.rows() - 1, nodes.cols());
        t_all.conservativeResize(t_all.rows() + t_sub.rows(), NoChange);
        t_all.bottomRows(t_sub.rows()) = t_sub;

        Y_all.push_back(Y);
        D_all.push_back(D);
        v0 = y.row(y.rows() - 1).transpose();
        nodes += VectorXd::Ones(t_step.size()) * dt;
    }
    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "jnfk_uniform run time: " << duration.count() << " s" << endl;
    }
    return {t_all, y_all, Y_all, D_all};
}

tuple<VectorXd, MatrixXd, vector<vector<MatrixXd>>, vector<vector<MatrixXd>>, vector<double>>
jfnk_adaptive(
    function<VectorXd(double, const VectorXd &)> f_eval,
    double t_init,
    double t_final,
    double dt_init,
    const Params &p,
    const VectorXd &y0,
    double tol,
    int n_iter_max_newton,
    double be_tol,
    double sdc_tol,
    bool do_print)
{
    auto start = std::chrono::high_resolution_clock::now();
    assert(t_final > t_init);
    assert(dt_init <= t_final - t_init);

    double h = dt_init;
    int n_nodes = p.n_nodes;
    int node_type = p.node_type;

    MatrixXd y_all;                 // 所有节点的解，维度[t_init~t_final所有节点数(每个时间步Δt中的n_nodes总和) * 维度m]
    VectorXd t_all;                 // 所有节点的时间，维度（t_init~t_final所有节点数）
    vector<vector<MatrixXd>> Y_all; // 每个时间步的，每次迭代的解，维度[时间步Δt数 * 迭代次数 * Δt内节点数n_nodes * 维度m]
    vector<vector<MatrixXd>> D_all; // 每个时间步的，每次迭代的误差，维度[时间步Δt数 * 迭代次数 * n_nodes * 维度m]
    vector<double> h_all;           // 每个时间步的步长，维度（t_init~t_final所有时间步Δt数）

    y_all.conservativeResize(1, y0.size());
    y_all.row(0) = y0.transpose();
    t_all.conservativeResize(1, 1);
    t_all(0) = t_init;

    double t0 = t_init;
    VectorXd v0 = y0;
    bool do_stop = false;

    double scaler_big_enough = 1.5;
    double scaler;
    bool do_return;

    int j = 0;
    // cout << "adaptive run" << endl;
    while (t0 < t_final && !do_stop)
    {
        assert(j < n_iter_max_newton);
        cout << "j: " << j << endl;

        VectorXd t1, t2;
        MatrixXd y1, y2;
        vector<vector<MatrixXd>> Y1, Y2;
        vector<vector<MatrixXd>> D1, D2;
        // cout << "t0: " << t0 << endl;
        // cout << "h: " << h << endl;
        tie(t1, y1, Y1, D1) = jfnk_uniform(f_eval, t0, t0 + h, 1, p, v0, n_iter_max_newton, be_tol, sdc_tol);
        // cout << "-------t1信息：---------- " << endl;
        // cout << t1 << endl;
        // cout << "-------y1信息：---------- " << endl;
        // cout << y1 << endl;
        tie(t2, y2, Y2, D2) = jfnk_uniform(f_eval, t0, t0 + h, 2, p, v0, n_iter_max_newton, be_tol, sdc_tol);
        // cout << "t2信息： " << endl;
        // cout << t2 << endl;
        // cout << "y2信息： " << endl;
        // cout << y2 << endl;

        scaler = step_size_scaler(y1.row(y1.rows() - 1).transpose(), y2.row(y2.rows() - 1).transpose(), n_nodes, 2, tol, node_type);
        // cout << "scaler: " << scaler << endl;
        scaler = adjust_scaler(scaler, scaler_big_enough);
        // cout << "scaler: " << scaler << endl;

        do_return = ((scaler < 1) || (scaler > scaler_big_enough));
        // cout << "do_return: " << do_return << endl;

        VectorXd t;
        MatrixXd y;
        MatrixXd y_sub;
        VectorXd t_sub;
        vector<vector<MatrixXd>> Y;
        vector<vector<MatrixXd>> D;

        if (do_return)
        {
            h = update_step_size(scaler * h, t0, t_final);
            // cout << "--------do_return------" << endl;
            // cout << "do_return h: " << h << endl;
            // cout << "do_return t0: " << t0 << endl;
            // cout << "do_return t0 + h: " << t0 + h << endl;
            // cout << "do_return v0: " << endl;
            // cout << v0 << endl;
            tie(t, y, Y, D) = jfnk_uniform(f_eval, t0, t0 + h, 1, p, v0, n_iter_max_newton, be_tol, sdc_tol);
            // cout << "do_return t信息： " << endl;
            // cout << t << endl;
            // cout << "do_return y信息： " << endl;
            // cout << y << endl;
        }
        else
        {
            t = t1;
            y = y1;
            Y = Y1;
            D = D1;
        }

        // cout << "t信息： " << endl;
        // cout << t << endl;
        // cout << "y信息： " << endl;
        // cout << y << endl;

        if (p.t(t.size() - 1) == 1.0)
        {

            y_sub = y.block(1, 0, y.rows() - 1, y.cols());
            y_all.conservativeResize(y_all.rows() + y_sub.rows(), NoChange);
            y_all.bottomRows(y_sub.rows()) = y_sub;

            t_sub = t.block(1, 0, t.rows() - 1, t.cols());
            t_all.conservativeResize(t_all.rows() + t_sub.rows(), NoChange);
            t_all.bottomRows(t_sub.rows()) = t_sub;
        }
        else
        {
            y_sub = y.block(1, 0, y.rows(), y.cols());
            y_all.conservativeResize(y_all.rows() + y_sub.rows(), NoChange);
            y_all.bottomRows(y_sub.rows()) = y_sub;

            t_sub = t.block(1, 0, t.rows(), t.cols());
            t_all.conservativeResize(t_all.rows() + t_sub.rows(), NoChange);
            t_all.bottomRows(t_sub.rows()) = t_sub;
        }

        if (j==100) {
            cout << "j: " << j << endl;
            cout << t0 << " / " << t_final << endl;
            auto cur = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = cur - start;
            cout << "adaptive run time: " << duration.count() << " s" << endl;
            cout << "y信息： " << endl;
            cout << y << endl;
            cout << "y_all信息： " << endl;
            cout << y_all << endl;
            cout << endl;
        }

        h_all.push_back(h);
        Y_all.push_back(Y[0]);
        D_all.push_back(D[0]);
        v0 = y.row(y.rows() - 1).transpose();
        t0 += h;

        do_stop = stopping_criteria(t_final, t0);
        // cout << "do_stop: " << do_stop << endl;
        h = update_step_size(h, t0, t_final);
        // cout << "h: " << h << endl;
        j++;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    cout << "adaptive run time: " << duration.count() << " s" << endl;
    if (do_print)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        cout << "adaptive run time: " << duration.count() << " s" << endl;
    }
    return {t_all, y_all, Y_all, D_all, h_all};
}
