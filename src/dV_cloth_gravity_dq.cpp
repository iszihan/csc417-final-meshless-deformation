#include <dV_cloth_gravity_dq.h>
#include <iostream>
void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    Eigen::VectorXd stacked_g(M.rows());
    for(int i=0; i<M.rows()/3; ++i){
        stacked_g(i*3) = g(0);
        stacked_g(i*3+1) = g(1);
        stacked_g(i*3+2) = g(2);
    }
    fg = - M * stacked_g;
}
