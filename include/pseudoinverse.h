#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>

// method for calculating the pseudo-Inverse as recommended by Eigen developers
Eigen::MatrixXd pseudoinverse(const Eigen::MatrixXd &a, double epsilon);