#include "preprocess.h"
#include <Eigen/Eigenvalues>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <stdexcept>
#include "points.h"

#define GAUSS_LEGENDRE 1
#define GAUSS_LOBATTO 2
#define GAUSS_RADAU 3
#define GAUSS_RADAU_2A 4

using namespace std;
using namespace boost::math::quadrature;

// 阶乘
long long factorial(int n)
{
    long long res = 1;
    for (int i = 2; i <= n; ++i)
        res *= i;
    return res;
}

// 数值积分（使用 Boost Gauss-Kronrod 算法）
IntegralResult integrate_quad(const std::function<double(double)> &f, double a, double b)
{

    gauss_kronrod<double, 21> gk;

    double error_estimate;
    // 积分计算：f为被积函数，[a,b]为积分区间，1e-12为绝对误差容忍度，1000为最大迭代次数
    double result = gk.integrate(f, a, b, 1e-12, 1000, &error_estimate);

    return {result, error_estimate}; // 返回积分结果和误差估计
}

// // 数值积分（简化版，实际应用建议用 Boost 或 Gauss-Kronrod 算法）
// IntegralResult integrate_quad(const std::function<double(double)>& f, double a, double b) {
//     // 此处用自适应辛普森积分示例（精度较低，需根据需求优化）
//     const int max_iter = 1000;
//     const double tol = 1e-8;
//     double h = b - a;
//     double c = (a + b) / 2;
//     double fa = f(a), fb = f(b), fc = f(c);
//     double S = (h/6) * (fa + 4*fc + fb);
//     double error = 0;

//     for (int i = 0; i < max_iter; ++i) {
//         double d = (a + c)/2, e = (c + b)/2;
//         double fd = f(d), fe = f(e);
//         double Sleft = (h/12) * (fa + 4*fd + 2*fc + 4*fe + fb);
//         double Sright = (h/12) * (fa + 4*fd + 2*fc + 4*fe + fb); // 实际应分左右区间计算
//         error = std::abs(Sleft + Sright - S) / 15;
//         S = (Sleft + Sright + S)/15;
//         h /= 2;
//         if (error < tol) break;
//     }
//     return {S, error};
// }

// 向后欧拉积分矩阵
Eigen::MatrixXd backward_euler_matrix(const Eigen::VectorXd &t)
{
    bool do_end_point;
    if (t[t.size() - 1] == 1.0)
    {
        do_end_point = true;
    }
    else
    {
        do_end_point = false;
    }

    Eigen::VectorXd tau;
    if (!do_end_point)
    {
        tau.resize(t.size() + 1);
        tau.head(t.size()) = t;
        tau(t.size()) = 1.0;
    }
    else
    {
        tau = t;
    }

    int n = tau.size();
    Eigen::MatrixXd S_p(n, n);
    S_p.setZero();

    /*
    S_p = {t0,
           t0,t1-t0,
           t0,t1-t0,t2-t1
           t0,t1-t0,t2-t1,t3-t2
           t0,t1-t0,t2-t1,t3-t2,t4-t3}
    */

    for (int c = 0; c < n; ++c)
    {
        if (c == 0)
        {
            S_p.col(c).setConstant(tau(c));
        }
        else
        {
            double delta = tau(c) - tau(c - 1);
            S_p.block(c, c, n - c, 1).setConstant(delta); // 下三角部分赋值
        }
    }
    return S_p;
}

// Gauss-Lobatto 节点生成（需补充具体实现，此处为占位符）
Eigen::VectorXd eigen_gauss_lobatto(int n)
{
    // 实际实现需计算 Legendre 多项式零点，此处简化返回等间隔节点（示例）
    vector<double> nodes_vec = gauss_lobatto(n);
    Eigen::VectorXd nodes(nodes_vec.size());
    for (int i = 0; i < nodes_vec.size(); ++i)
    {
        nodes(i) = nodes_vec[i];
    }
    return nodes;
}

// Gauss-Lobatto 误差常数,预估数值积分的误差范围
double gauss_lobatto_error_constant(int n)
{
    double top = n * pow(n - 1, 3) * pow(factorial(n - 2), 4);
    double bot = (2 * n - 1) * pow(factorial(2 * n - 2), 3);
    return top / bot;
}

// 生成指定类型节点（需补充其他节点类型实现）
Eigen::VectorXd get_nodes(int n, int node_type)
{
    vector<double> nodes_vec;
    Eigen::VectorXd nodes;
    // 此处仅实现 Gauss-Lobatto 示例，其他节点类型需补充
    switch (node_type)
    {
    case GAUSS_LOBATTO:
        nodes_vec = gauss_lobatto(n);
        nodes = Eigen::Map<Eigen::VectorXd>(nodes_vec.data(), nodes_vec.size());
        break;
    case GAUSS_LEGENDRE:
        nodes_vec = gauss_legendre(n);
        nodes = Eigen::Map<Eigen::VectorXd>(nodes_vec.data(), nodes_vec.size());
        break;
    case GAUSS_RADAU:
        nodes_vec = gauss_radau(n);
        nodes = Eigen::Map<Eigen::VectorXd>(nodes_vec.data(), nodes_vec.size());
        break;
    case GAUSS_RADAU_2A:
        nodes_vec = gauss_radau_2a(n);
        nodes = Eigen::Map<Eigen::VectorXd>(nodes_vec.data(), nodes_vec.size());
        break;
    default:
        cout << "get_node_error" << endl;
        break;
    }

    // 转换到 [0, 1] 区间（若原节点在 [-1, 1]）
    nodes = (nodes.array() + 1) / 2.0;
    return nodes;
}

// 谱积分矩阵 S
Eigen::MatrixXd spectral_matrix(const Eigen::VectorXd &t)
{
    int n = t.size();
    bool do_end_point;
    int rows;
    if (t(n - 1) == 1.0)
    {
        do_end_point = true;
    }
    else
    {
        do_end_point = false;
    }

    if (do_end_point)
    {
        rows = n;
    }
    else
    {
        rows = n + 1;
    }

    Eigen::MatrixXd S(rows, n);
    S.setZero();

    for (int c = 0; c < n; c++)
    { // 列索引：拉格朗日基函数序号
        // 构造拉格朗日基函数分子：f(x) = product_{k≠c} (x - t[k])
        auto f = [&](double x)
        {
            double prod = 1.0;
            for (int k = 0; k < n; ++k)
            {
                if (k != c)
                {
                    prod *= (x - t(k));
                }
            }
            return prod;
        };
        // 分母：denom = product_{k≠c} (t[c] - t[k])
        double denom = 1.0;
        for (int k = 0; k < n; ++k)
        {
            if (k != c)
            {
                denom *= (t(c) - t(k));
            }
        }

        // 计算积分：S[r][c] = ∫₀^{t[r]} f(x)dx / denom
        for (int r = 0; r < n; ++r)
        { // 行索引：积分上限 t[r]
            IntegralResult res = integrate_quad(f, 0, t(r));
            S(r, c) = res.value / denom;
        }
        // 若不包含终点，额外计算积分上限为 1 的行
        if (!do_end_point)
        {
            IntegralResult res = integrate_quad(f, 0, 1.0);
            S(rows - 1, c) = res.value / denom;
        }
    }
    return S;
}

// 谱半径计算
double spectral_radius(int node_type, const Eigen::MatrixXd &S, const Eigen::MatrixXd &S_p)
{
    Eigen::MatrixXd C;
    int n = S.rows();

    switch (node_type)
    {
    case GAUSS_LEGENDRE:
    {
        int size = n - 1;
        Eigen::MatrixXd S_p_sub = S_p.block(0, 0, size, size); // S_p[1:, 1:]
        Eigen::MatrixXd S_sub = S.block(0, 0, size, size);     // S[1:, 1:]
        C = Eigen::MatrixXd::Identity(size, size) - S_p_sub.inverse() * S_sub;
        break;
    }
    case GAUSS_LOBATTO:
    {
        int size = n - 1;
        Eigen::MatrixXd S_p_sub = S_p.block(1, 1, size, size); // S_p[1:, 1:]
        Eigen::MatrixXd S_sub = S.block(1, 1, size, size);     // S[1:, 1:]
        C = Eigen::MatrixXd::Identity(size, size) - S_p_sub.inverse() * S_sub;
        break;
    }
    case GAUSS_RADAU:
    {
        int size = n - 1;
        Eigen::MatrixXd S_p_sub = S_p.block(1, 1, size, size); // S_p[1:, 1:]
        Eigen::MatrixXd S_sub = S.block(1, 0, size, size);     // S[1:, 1:]
        C = Eigen::MatrixXd::Identity(size, size) - S_p_sub.inverse() * S_sub;
        break;
    }
    case GAUSS_RADAU_2A:
    {
        int size = n;
        C = Eigen::MatrixXd::Identity(size, size) - S_p.inverse() * S;
        break;
    }
    // 其他节点类型（GAUSS_LEGENDRE 等）需补充实现
    default:
        cout << "spectral_radius_error" << endl;
        break;
    }

    Eigen::EigenSolver<Eigen::MatrixXd> es(C);
    return es.eigenvalues().cwiseAbs().maxCoeff();
}