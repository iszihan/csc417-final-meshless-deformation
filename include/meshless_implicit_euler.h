#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>
#include <igl/polar_dec.h>
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE> 
inline void meshless_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,
                                    Eigen::MatrixXd &V0, const Eigen::Vector3d &center_of_mass0,
                                    FORCE &force, Eigen::VectorXd &tmp_force) {
    std::cout<<"inside meshless update function"<<std::endl;

    //gather forces
    force(tmp_force,q,qdot);
    
    //update without the goal position fitting
    Eigen::VectorXd qdot_tmp = qdot + dt * tmp_force/mass;     
    Eigen::VectorXd q_tmp = q + dt * qdot_tmp;
    //Eigen::VectorXd q2_tmp = q + dt * qdot_tmp;
    // q2_tmp(18) = q(18);
    // q2_tmp(19) = q(19);
    // q2_tmp(20) = q(20);
    // q_tmp = q2_tmp;

    //compute goal positions
    //get current center of mass
    Eigen::Vector3d center_of_masst;
    Eigen::MatrixXd Vt = Eigen::Map<Eigen::MatrixXd>(q_tmp.data(),3,q_tmp.rows()/3);
    center_of_masst = Vt.transpose().colwise().mean();
    
    // std::cout<<center_of_masst<<std::endl;
    // if(isnan(center_of_masst(0))){
    //     std::exit(1);
    // }

    //get q: vertex position relative to CoM at t0
    Eigen::MatrixXd P = V0.rowwise() - center_of_mass0.transpose();
    //get p: vertex position relative to current CoM 
    Eigen::MatrixXd Q = Vt.transpose().rowwise() - center_of_masst.transpose();
    Eigen::Matrix3d Apq = mass * P.transpose() * Q; // 3 x 3 
    Eigen::Matrix3d Aqq = mass * Q.transpose() * Q; // 3 x 3
    
    //polar decomposition to get rotation
    
    // 1. formula in paper
    // Eigen::Matrix3d S = Apq.transpose() * Apq;
    // S = S.cwiseSqrt(); //this might have some assumption on A?
    // Eigen::Matrix3d R = Apq * S.inverse();

    // 2. igl method that works
    Eigen::Matrix3d R;
    Eigen::Matrix3d S;
    igl::polar_dec(Apq, R, S);
    
    // 3. use svd
    // w, s, vh = svd(a, full_matrices=False)
    // u = w.dot(vh)
    // if side == 'right':
    //     # a = up
    //     p = (vh.T.conj() * s).dot(vh)
    // else:
    //     # a = pu
    //     p = (w * s).dot(w.T.conj())
    // return u, p
    //Eigen::JacobiSVD<MatrixXd> svd(Apq, ComputeThinU | ComputeThinV);
    
    Eigen::MatrixXd rotatedP = (R * P.transpose()).transpose();
    Eigen::MatrixXd gt = rotatedP.rowwise() + center_of_masst.transpose();
    // std::cout<<"S:"<<std::endl;
    // std::cout<<S<<std::endl;
    // std::cout<<"R:"<<std::endl;
    // std::cout<<R<<std::endl;
    std::cout<<"qdot_tmp:"<<std::endl;
    std::cout<<qdot_tmp<<std::endl;
    std::cout<<"q_tmp:"<<std::endl;
    std::cout<<q_tmp<<std::endl;
    std::cout<<"goals:"<<std::endl;
    std::cout<<gt<<std::endl;
    
    //flatten gt
    Eigen::VectorXd gt_flatten;
    gt_flatten.resize(gt.rows()*gt.cols());
    Eigen::MatrixXd gtt = gt.transpose();
    gt_flatten = Eigen::Map<Eigen::VectorXd>(gtt.data(), gtt.rows()*gtt.cols());
    
    //update
    double alpha = 1.0;
    qdot = qdot_tmp + alpha * ((gt_flatten-q_tmp))/dt;     
    q = q_tmp + dt * qdot;
    //Eigen::VectorXd q = q_tmp + dt * qdot;
    // q2(18) = q(18);
    // q2(19) = q(19);
    // q2(20) = q(20);
    // q = q2;
}
