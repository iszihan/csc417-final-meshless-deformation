#include <rigid_body_jacobian.h>

void rigid_body_jacobian(Eigen::Matrix36d &J, 
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, 
                         Eigen::Ref<const Eigen::Vector3d> x) {

    // R [[X]^T I] [R^t  0 ] = [R[X]^TR^T I]
    //             [ 0  R^t]  
    Eigen::Matrix3d X_cp;
    X_cp<<0,-x(2),x(1),x(2),0,-x(0),-x(1),x(0),0;
    Eigen::Matrix36d XI;
    XI.block(0,0,3,3)=X_cp.transpose();
    XI.block(0,3,3,3)=Eigen::Matrix3d::Identity();
    Eigen::Matrix66d A;
    A.setZero(6,6);
    A.block(0,0,3,3) = R.transpose();
    A.block(3,3,3,3) = R.transpose();
    J = R * XI * A;
    
}