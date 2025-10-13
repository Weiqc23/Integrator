#ifndef POINTS_H
#define POINTS_H

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/tools/roots.hpp>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <set>

// 节点类型常量定义
#define GAUSS_LEGENDRE 1
#define GAUSS_LOBATTO 2
#define GAUSS_RADAU 3
#define GAUSS_RADAU_2A 4

// 求解多项式根（输入：多项式函数；输出：根的有序列表）
std::vector<double> find_roots(const std::function<double(double)>& poly);

// 生成 Gauss-Legendre 节点（n 个节点）
std::vector<double> gauss_legendre(int n);

// 生成 Gauss-Lobatto 节点（n 个节点，包含端点 ±1）
std::vector<double> gauss_lobatto(int n);

std::vector<double> gauss_radau(int n);

std::vector<double> gauss_radau_2a(int n);

// 计算 Legendre 多项式 P_n(x) 的导数 P'_n(x)
double legendre_deriv(int n, double x);

#endif // POINTS_H