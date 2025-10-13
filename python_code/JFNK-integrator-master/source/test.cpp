#include <Eigen/Dense>
#include "points.h"
#include <iostream>

int main() {
    int n = 5;
    // 步骤1：先获取 std::vector<double> 类型的节点
    std::vector<double> nodes_vec = gauss_lobatto(n);
    // 步骤2：通过 Eigen::Map 将 std::vector<double> 转换为 Eigen::VectorXd
    Eigen::VectorXd nodes_cpp = Eigen::Map<Eigen::VectorXd>(nodes_vec.data(), nodes_vec.size());
    
    std::cout << "C++ nodes: " << nodes_cpp.transpose() << std::endl;
    return 0;
}