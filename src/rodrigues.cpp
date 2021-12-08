#include <rodrigues.h>
#include <cmath>

void rodrigues(Eigen::Matrix3d &R, Eigen::Ref<const Eigen::Vector3d> omega) {

    Eigen::Vector3d axis = omega.normalized();
    double theta = omega.norm();
    Eigen::Matrix3d axis_cp;
    axis_cp<<0,-axis(2),axis(1),axis(2),0,-axis(0),-axis(1),axis(0),0;
    R = Eigen::Matrix3d::Identity() + std::sin(theta)*axis_cp + (1 - std::cos(theta))*axis_cp*axis_cp;

}