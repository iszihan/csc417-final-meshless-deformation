#include <dV_spring_particle_particle_dq.h>
#include <iostream>
void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    //std::cout<<"stiffness"<<std::endl;
    //std::cout<<stiffness<<std::endl;
    //std::cout<<(l0 - (q0-q1).norm())<<std::endl;
    //std::cout<<(q0-q1)<<std::endl;   
    Eigen::Vector3d f_q0 = - stiffness * (q0-q1) * (l0 - (q0-q1).norm()) / (q0-q1).norm(); 
    Eigen::Vector3d f_q1 =  stiffness * (q0-q1) * (l0 - (q0-q1).norm()) / (q0-q1).norm();
    f(0) = f_q0(0);
    f(1) = f_q0(1);
    f(2) = f_q0(2);
    f(3) = f_q1(0);
    f(4) = f_q1(1);
    f(5) = f_q1(2);
}