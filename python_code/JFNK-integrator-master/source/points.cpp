#include "points.h"
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/tools/roots.hpp>
#include <vector>
#include <functional>
#include <set>

using namespace std;
using namespace boost::math;
using namespace boost::math::tools;

#define GAUSS_LEGENDRE 1
#define GAUSS_LOBATTO 2
#define GAUSS_RADAU 3
#define GAUSS_RADAU_2A 4

// 输入：多项式poly，多项式次数degree
// 输出：多项式poly=0的根
vector<double> find_roots(const function<double(double)> &poly)
{
    set<double> roots_set;

    // --------------------------
    // 步骤1：生成初始根猜测（采样法）
    // --------------------------
    vector<double> initial_guesses;
    const double step = 0.01;
    double x_pre = -1.0;
    double f_pre = poly(x_pre);

    // 检测根的大致位置
    for (double x = x_pre + step; x <= 1.0 + step; x += step)
    {
        double f_cur = poly(x);

        if (f_pre * f_cur <= 0)
        {
            if (fabs(f_cur) < 1e-12)
            {
                initial_guesses.push_back(x);
            }
            else if (fabs(f_pre) < 1e-12)
            {
                initial_guesses.push_back(x_pre);
            }
            else
            {
                double x0 = x_pre - f_pre * (x - x_pre) / (f_cur - f_pre);
                initial_guesses.push_back(x0);
            }
        }

        x_pre = x;
        f_pre = f_cur;
    }

    for (double x0 : initial_guesses)
    {
        auto f = [&](double z)
        { return poly(z); };

        double a = x0 - 0.1; // 初始左边界
        double b = x0 + 0.1; // 初始右边界
        double fa = f(a);
        double fb = f(b);

        int expand_count = 0;
        while (fa * fb > 0 && expand_count < 10)
        { // 最多扩展10次
            a -= 0.1;
            b += 0.1;
            fa = f(a);
            fb = f(b);
            expand_count++;
        }

        if (fa * fb > 0)
        { // 仍无符号变化，跳过该猜测
            continue;
        }

        // 使用Boost的TOMS748算法求解根
        const int precision_digits = 50;
        boost::math::tools::eps_tolerance<double> tol(precision_digits); // 创建容差对象
        uintmax_t max_iter = 1000;                                       // 设置最大迭代次数

        // 调用toms748_solve求解根
        pair<double, int> root = boost::math::tools::toms748_solve(f, a, b, fa, fb, tol, max_iter);

        if (root.second == 0)
        {
            roots_set.insert(root.first);
        }
    }

    return vector<double>(roots_set.begin(), roots_set.end());
}

// 计算Legendre多项式P_n(x)-->boost库
// double legendre_poly(int n, double x) {
//     return boost::math::legendre_p(n, x);  // Boost.Math的Legendre多项式实现
// }

// 计算Legendre多项式的导数P'_n(x)
double legendre_deriv(int n, double x)
{
    if (n == 0)
        return 0.0; // P_0(x)=1，导数为0

    if (x < -1.0)
        x = -1.0;
    else if (x > 1.0)
        x = 1.0;
    // 利用递推公式：P'_n(x) = [n*x*P_n(x) - n*P_{n-1}(x)] / (x² - 1)
    double p_n = boost::math::legendre_p(n, x);
    double p_n_1 = boost::math::legendre_p(n - 1, x);
    return (n * x * p_n - n * p_n_1) / (x * x - 1 + 1e-12); // 避免除零
}

// Gauss-Legendre节点：n次Legendre多项式P_n(x)的根
vector<double> gauss_legendre(int n)
{
    if (n <= 0)
        return {};

    std::function<double(double)> f = [n](double x)
    { return boost::math::legendre_p(n, x); };
    return find_roots(f); // n次多项式最多n个根
}

// Gauss-Legendre节点：n-1次Legendre多项式P_n(x)的导数
vector<double> gauss_lobatto(int n)
{
    if (n <= 1)
        return {-1.0, 1.0}; // 至少2个节点（端点）
    // 定义函数f(x) = P'_{n-1}(x)
    auto f = [n](double x)
    { return legendre_deriv(n - 1, x); };
    vector<double> roots = find_roots(f); // P'_{n-1}是n-2次多项式
    // 添加端点并去重排序
    roots.push_back(-1.0);
    roots.push_back(1.0);
    sort(roots.begin(), roots.end());
    return roots;
}

// Gauss-Radau节点：P_n(x) + P_{n-1}(x)的根（包含左端点-1）
vector<double> gauss_radau(int n)
{
    if (n <= 0)
        return {};
    // 定义函数f(x) = P_n(x) + P_{n-1}(x)
    auto f = [n](double x)
    {
        return boost::math::legendre_p(n, x) + legendre_p(n - 1, x);
    };
    vector<double> roots = find_roots(f);
    sort(roots.begin(), roots.end());
    return roots;
}

// Gauss-Radau 2A节点：P_n(x) - P_{n-1}(x)的根（包含右端点1）
vector<double> gauss_radau_2a(int n)
{
    if (n <= 0)
        return {};
    // 定义函数f(x) = P_n(x) - P_{n-1}(x)
    auto f = [n](double x)
    {
        return boost::math::legendre_p(n, x) - legendre_p(n - 1, x);
    };
    vector<double> roots = find_roots(f);
    sort(roots.begin(), roots.end());
    return roots;
}