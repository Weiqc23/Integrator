#ifndef PREPROCESS_H
#define PREPROCESS_H

#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#define GAUSS_LEGENDRE 1
#define GAUSS_LOBATTO 2
#define GAUSS_RADAU 3
#define GAUSS_RADAU_2A 4

using namespace boost::math::quadrature;
using namespace std;

struct IntegralResult {
    double value;
    double error;
    IntegralResult(double v, double e) : value(v), error(e) {}
};

long long factorial(int n);

IntegralResult integrate_quad(const std::function<double(double)> &f, double a, double b);

Eigen::MatrixXd backward_euler_matrix(const Eigen::VectorXd &t);

Eigen::VectorXd eigen_gauss_lobatto(int n);

double gauss_lobatto_error_constant(int n);

Eigen::VectorXd get_nodes(int n, int node_type);

Eigen::MatrixXd spectral_matrix(const Eigen::VectorXd &t);

double spectral_radius(int node_type, const Eigen::MatrixXd &S, const Eigen::MatrixXd &S_p);

#endif // PREPROCESS_H