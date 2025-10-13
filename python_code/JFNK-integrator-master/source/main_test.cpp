#include "integrator_liner.h"
#include "points.h"
#include "matplotlibcpp.h"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <chrono>

namespace plt = matplotlibcpp;
using namespace std;
using namespace Eigen;

// VectorXd f_eval(double t, const VectorXd &y)
// {
//     return -y;
// }

VectorXd loadCSVVectorXd(const string &path)
{
    ifstream file(path);
    string line;
    vector<double> vec;
    int rows;

    rows = 0;
    while (getline(file, line))
    {
        vec.push_back(stod(line));
        rows++;
    }

    VectorXd vecd(rows);
    for (int i = 0; i < rows; i++)
    {
        vecd(i) = vec[i];
    }

    return vecd;
}

MatrixXd loadCSVMatrix(const string &path)
{
    ifstream file(path);
    string line;
    vector<vector<double>> mat_vec;
    int rows, cols;

    rows = 0;
    while (getline(file, line))
    {
        // cout << "line: " << line << endl;
        stringstream ss(line);
        string cell;
        vector<double> row;
        while (getline(ss, cell, ','))
        {
            // cout << "cell: " << cell << endl;
            row.push_back(stod(cell));
        }
        mat_vec.push_back(row);
        cols = mat_vec[rows].size();
        rows++;
    }

    MatrixXd mat(rows, cols);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            mat(i, j) = mat_vec[i][j];
        }
    }

    return mat;
}

MatrixXd u_func(double t)
{
    MatrixXd u(1, 1);

    double peak_amplitude = 735e3 / sqrt(3.0) * sqrt(2.0);
    double freq = 60.0; // Hz
    double w = 2.0 * M_PI * freq;
    double phase = 0.0; // rad

    u(0, 0) = peak_amplitude * sin(w * t + phase);

    return u;
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    double t_init = 0.0;
    double t_final = 0.1;
    double dt_init = 100e-6;
    double tol = 1e-6;

    MatrixXd A = loadCSVMatrix("C:\\Users\\wei13\\Desktop\\Basic_Code\\Basic_Code\\Ex_AC_System\\problem_dependent\\A.csv");
    MatrixXd B = loadCSVMatrix("C:\\Users\\wei13\\Desktop\\Basic_Code\\Basic_Code\\Ex_AC_System\\problem_dependent\\B.csv");
    VectorXd y0 = loadCSVVectorXd("C:\\Users\\wei13\\Desktop\\Basic_Code\\Basic_Code\\Ex_AC_System\\problem_dependent\\y0.csv");

    // cout << "A: " << endl << A << endl;
    // cout << "B: " << endl << B << endl;
    // cout << "y0: " << endl << y0 << endl;
    auto u = [&](double t)
    {
        return u_func(t);
    };

    auto f_eval = [&](double t, const VectorXd &y)
    {
        VectorXd res = A * y + B * u_func(t);
        return res;
    };
    // 构造参数对象（3 个节点，高斯-Lobatto）
    
    Params p(4, y0.size(), GAUSS_LOBATTO);
    // cout << "params.n_nodes: " << p.n_nodes << endl;
    // cout << "params.node_type: " << p.node_type << endl;
    // cout << "params.t: " << p.t << endl;
    // 调用 C++ 版 jfnk_adaptive
    // auto [t_all, y_all, Y_all, D_all, h_all] = jfnk_adaptive(
    //     f_eval, t_init, t_final, dt_init, p, y0, tol);

    auto [t_all, y_all, Y_all, D_all, h_all] = jfnk_adaptive(
        f_eval, A, B, u, t_init, t_final, dt_init, p, y0, tol);

    // 打印t_all信息
    int t_size = t_all.size();
    cout << "\nt_all信息：" << endl;
    cout << "t_all维度: " << t_size << endl;
    // cout << "t_all内容: " << endl
    //      << t_all.transpose() << endl; // 转置为行向量便于查看

    // 打印y_all信息
    double y_rows = y_all.rows();
    double y_cols = y_all.cols();
    cout << "\ny_all信息：" << endl;
    cout << "y_all维度: " << y_rows << " x " << y_cols << endl;
    cout << "y_all.row(0)内容: " << endl;
    cout << fixed << setprecision(10) << y_all.row(0) << endl;
    // cout << y_all.row(0) << endl;
    cout << "y_all.row(-1)内容: " << endl;
    cout << fixed << setprecision(10) << y_all.row(y_all.rows() - 1) << endl;

    vector<double> t_vec(t_all.data(), t_all.data() + t_all.size());
    vector<double> y_vec(y_all.rows());
    for (int i = 0; i < y_all.rows(); i++)
    {
        y_vec[i] = y_all(i, 0);
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "运行时间: " << duration.count() << " 秒" << endl;
    // 绘制图形

    plt::figure_size(1200, 780);
    plt::plot(t_vec, y_vec, "b-");
    plt::title("JFNK_cpp");
    plt::xlabel("t_all)");
    plt::ylabel("y_all(row0)");
    plt::grid(true);
    plt::show();

    return 0;
}
