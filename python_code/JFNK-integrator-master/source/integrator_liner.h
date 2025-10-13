#ifndef INTEGRATOR_LINER_H
#define INTEGRATOR_LINER_H

#include <Eigen/Dense>
#include <vector>
#include <functional>
#include <tuple>

// ---------------- 常量定义 ----------------
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

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

struct Params
{
    int n_nodes;            // 每个时间步的节点数
    int m;                  // 方程维度
    int node_type;          // 节点类型
    VectorXd t;             // 节点
    MatrixXd S;             // 谱矩阵
    MatrixXd S_p;           // 谱矩阵的导数
    double spectral_radius; // 谱半径

    Params(int n_nodes, int m, int node_type = GAUSS_LOBATTO);
};

// ---------------- 函数声明 ----------------

// 牛顿迭代求解器
VectorXd root(std::function<VectorXd(const VectorXd &)> f,
              const VectorXd &x0,
              double tol = BE_TOL,
              int max_iter = N_ITER_MAX_NEWTON,
              bool do_print = false);

// backward euler 单步
VectorXd backward_euler_node(function<VectorXd(double, const VectorXd &)> f_eval,
                             const MatrixXd &A,
                             const MatrixXd &B,
                             function<MatrixXd(double)> u,
                             double t,
                             double h,
                             const VectorXd &rhs,
                             const VectorXd &x0,
                             double be_tol = BE_TOL,
                             bool do_print = false);

// backward euler 整个时间步
MatrixXd backward_euler(function<VectorXd(double, const VectorXd &)> f_eval,
                        const MatrixXd &A,
                        const MatrixXd &B,
                        function<MatrixXd(double)> u,
                        const VectorXd &t,
                        const VectorXd &y0,
                        const MatrixXd &S_p,
                        double be_tol = BE_TOL,
                        bool do_print = false);

pair<VectorXd, VectorXd> sdc_node(function<VectorXd(double, const VectorXd &)> f_eval,
                                  const MatrixXd &A,
                                  const MatrixXd &B,
                                  function<MatrixXd(double)> u,
                                  double t,
                                  double h,
                                  const VectorXd &rhs,
                                  const VectorXd &y_old,
                                  double be_tol = BE_TOL,
                                  bool do_print = false);

tuple<MatrixXd, MatrixXd>
sdc_sweep(function<VectorXd(double, const VectorXd &)> f_eval,
          const MatrixXd &A,
          const MatrixXd &B,
          function<MatrixXd(double)> u,
          const VectorXd &t,
          const VectorXd &y0,
          const MatrixXd &y_old,
          const MatrixXd &F,
          const MatrixXd &S,
          const MatrixXd &S_p,
          double be_tol = BE_TOL,
          bool do_print = false);

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool>
sdc(function<VectorXd(double, const VectorXd &)> f_eval,
    const MatrixXd &A,
    const MatrixXd &B,
    function<MatrixXd(double)> u,
    const VectorXd &t,
    const VectorXd &y0,
    const MatrixXd &y_old,
    const MatrixXd &S,
    const MatrixXd &S_p,
    int n_iter_max_sdc = N_ITER_MAX_SDC,
    double be_tol = BE_TOL,
    double sdc_tol = SDC_TOL,
    bool do_print = false);

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool, bool, VectorXd>
jfnk_initial(function<VectorXd(double, const VectorXd &)> f_eval,
             const MatrixXd &A,
             const MatrixXd &B,
             function<MatrixXd(double)> u,
             const VectorXd &t,
             const VectorXd &y0,
             const MatrixXd &y_approx,
             const MatrixXd &S,
             const MatrixXd &S_p,
             double spectral_radius,
             double be_tol = BE_TOL,
             double sdc_tol = SDC_TOL,
             int n_iter_max = N_ITER_MAX_SDC,
             bool do_print = false);

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool>
jfnk_iterations(function<VectorXd(double, const VectorXd &)> f_eval,
                const MatrixXd &A,
                const MatrixXd &B,
                function<MatrixXd(double)> u,
                const VectorXd &t,
                const VectorXd &y0,
                const MatrixXd &y_init,
                const MatrixXd &S,
                const MatrixXd &S_p,
                int n_iter_max_newton = N_ITER_MAX_NEWTON,
                double be_tol = BE_TOL,
                double sdc_tol = SDC_TOL,
                bool do_print = false);

tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool, bool, VectorXd>
jfnk(function<VectorXd(double, const VectorXd &)> f_eval,
     const MatrixXd &A,
     const MatrixXd &B,
     function<MatrixXd(double)> u,
     const VectorXd &t,
     const VectorXd &y0,
     const MatrixXd &y_approx,
     const MatrixXd &S,
     const MatrixXd &S_p,
     double spectral_radius,
     int n_iter_max_newton = N_ITER_MAX_NEWTON,
     double be_tol = BE_TOL,
     double sdc_tol = SDC_TOL,
     bool do_print = false);

// JFNK 单步
tuple<MatrixXd, vector<MatrixXd>, vector<MatrixXd>, bool, bool, VectorXd>
jfnk_step(function<VectorXd(double, const VectorXd &)> f_eval,
          const MatrixXd &A,
          const MatrixXd &B,
          function<MatrixXd(double)> u,
          const Params &p,
          double dt,
          const VectorXd &t,
          const VectorXd &y0,
          int n_iter_max_newton = N_ITER_MAX_NEWTON,
          double be_tol = BE_TOL,
          double sdc_tol = SDC_TOL,
          bool do_print = false);

// JFNK uniform 步长
tuple<VectorXd, MatrixXd, vector<vector<MatrixXd>>, vector<vector<MatrixXd>>>
jfnk_uniform(function<VectorXd(double, const VectorXd &)> f_eval,
             const MatrixXd &A,
             const MatrixXd &B,
             function<MatrixXd(double)> u,
             double t_init,
             double t_final,
             int n_steps,
             const Params &p,
             const VectorXd &y0,
             int n_iter_max_newton = N_ITER_MAX_NEWTON,
             double be_tol = BE_TOL,
             double sdc_tol = SDC_TOL,
             bool do_print = false);

// JFNK adaptive 步长
tuple<VectorXd, MatrixXd, vector<vector<MatrixXd>>, vector<vector<MatrixXd>>, vector<double>>
jfnk_adaptive(function<VectorXd(double, const VectorXd &)> f_eval,
              const MatrixXd &A,
              const MatrixXd &B,
              function<MatrixXd(double)> u,
              double t_init,
              double t_final,
              double dt_init,
              const Params &p,
              const VectorXd &y0,
              double tol,
              int n_steps_max = N_STEPS_MAX_ADAPTIVE,
              double be_tol = BE_TOL,
              double sdc_tol = SDC_TOL,
              bool do_print = false);

#endif // INTEGRATOR_H
