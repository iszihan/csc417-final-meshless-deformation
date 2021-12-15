#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>
#include <pseudoinverse.h>

// method for calculating the pseudo-Inverse as recommended by Eigen developers
 Eigen::MatrixXd pseudoinverse(const Eigen::MatrixXd &a, double epsilon)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(a ,Eigen::ComputeFullU | Eigen::ComputeFullV);
    // For a non-square matrix
    // Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}