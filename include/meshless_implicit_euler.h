#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>
#include <igl/polar_dec.h>
#include <pseudoinverse.h>
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
inline void 
meshless_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,
                                    Eigen::SparseMatrixd &P, Eigen::VectorXd &x0,
                                    std::vector<Eigen::Vector3d> centers_of_mass, std::vector<Eigen::MatrixXd> &Qs, 
                                    std::vector<std::vector<int>> clusters, Eigen::VectorXi vertex_cluster_count,
                                    int method, FORCE &force, Eigen::VectorXd &tmp_force, Eigen::Vector3d & comt, std::vector<Eigen::Matrix3d>& Sps) {

    //std::cout<<"inside integration..."<<std::endl;
    //gather forces
    force(tmp_force,q,qdot);
    //std::cout<<tmp_force.norm()<<std::endl;    

    //update all vertices without the goal position fitting
    Eigen::VectorXd qdot_tmp = qdot + dt * tmp_force/mass;     
    Eigen::VectorXd q_tmp = q + dt * qdot_tmp;
    //keep fixed points unchanged
    q_tmp = P.transpose() * P * q_tmp + x0;
    
    //std::cout<<"Meshless integration with cluster size = "<<clusters.size()<<std::endl;
    Eigen::VectorXd qdot_holder(qdot.rows());
    qdot_holder.setZero();
    for(int ic=0; ic<clusters.size(); ++ic){

        //compile qi qdoti for current cluster vertices
        Eigen::VectorXd qi_tmp;
        Eigen::VectorXd qdoti_tmp;
        qi_tmp.setZero(clusters.at(ic).size()*3);
        qdoti_tmp.setZero(clusters.at(ic).size()*3);
        for(int iv=0;iv<clusters.at(ic).size();++iv){
            qi_tmp.segment<3>(iv*3) = q_tmp.segment<3>(clusters.at(ic).at(iv)*3);
            qdoti_tmp.segment<3>(iv*3) = qdot_tmp.segment<3>(clusters.at(ic).at(iv)*3);
        }
        
        //compute goal positions
        //get current center of mass
        Eigen::Vector3d center_of_masst;
        comt = center_of_masst;
        Eigen::MatrixXd Vt = Eigen::Map<Eigen::MatrixXd>(qi_tmp.data(),3,qi_tmp.rows()/3);
        center_of_masst = Vt.transpose().colwise().mean();
        //get p and q: vertex position relative to current and original CoM 
        Eigen::MatrixXd _Q = Qs.at(ic);
        Eigen::MatrixXd _P = Vt.transpose().rowwise() - center_of_masst.transpose();
        Eigen::Matrix3d Aqq = (mass * _Q.transpose() * _Q).inverse(); // 3 x 3
        Eigen::Matrix3d Apq = mass * _P.transpose() * _Q; // 3 x 3 
        Eigen::Matrix3d A = Apq * Aqq;

        //polar decomposition to get rotation
        Eigen::Matrix3d R;
        Eigen::Matrix3d S;
        //std::cout<<"before poalr dec"<<std::endl;
        igl::polar_dec(Apq, R, S);

        //plasticity 
        if(method == 3){
            //plasticity related
            double c_yield = 0.002;
            double c_max = 0.008;
            double c_creep = 0.001;
            Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
            double deformation = (S - I).norm();
            if(deformation > c_yield){
                // std::cout<<"plasticity deformation:"<<deformation<<std::endl;
                Eigen::Matrix3d Sp = (I + c_creep * (S-I)) * Sps.at(ic);
                // std::cout<<"updated Sp:"<<Sp<<std::endl;
                // std::cout<<"Sp-I.norm:"<<(Sp-I).norm()<<std::endl;
                //bound by a max
                if((Sp-I).norm()>c_max){
                    // std::cout<<"bigger than maximum deformation"<<std::endl;
                    Sp = I + c_max * (Sp-I) / (Sp-I).norm();
                    // std::cout<<Sp<<std::endl;
                    // std::cout<<(Sp-I).norm()<<std::endl;
                }
                Sp = Sp / std::cbrt(Sp.determinant());
                // std::cout<<Sp<<std::endl;
                //update q 
                _Q = (Sp * _Q.transpose()).transpose();
                Aqq = (mass * _Q.transpose() * _Q).inverse(); // 3 x 3
                Apq = mass * _P.transpose() * _Q; // 3 x 3 
                A = Apq * Aqq;
                igl::polar_dec(Apq, R, S);
            }
        }
        //std::cout<<"end of plasticity stuff"<<std::endl;


        //get p: vertex position relative to current CoM 
        Eigen::MatrixXd transformedP;
        if(method == 0){
            //rigid
            //std::cout<<"using method 0: rigid"<<std::endl;
            Eigen::Matrix3d T;
            T = R;
            transformedP = (T * _Q.transpose()).transpose();
        }else if(method == 1 || method == 3){
            //linear
            //std::cout<<"using method 1: linear"<<std::endl;
            double beta = 0.7;
            Eigen::Matrix3d T;
            A = A / std::cbrt(A.determinant());
            T = beta * A + (1-beta) * R;
            transformedP = (T * _Q.transpose()).transpose();
        }else if(method == 2){
            //quadratic
            //std::cout<<"using method 2: quadratic"<<std::endl;
            
            //TODO: move this to pre-computation?
            Eigen::MatrixXd _Qdelta;
            _Qdelta.setZero(_Q.rows(),9);
            for(int r=0; r<_Q.rows(); ++r){
                Eigen::VectorXd _qdelta;
                _qdelta.setZero(9);
                _qdelta<<_Q.row(r)(0),
                            _Q.row(r)(1),
                            _Q.row(r)(2),
                            std::pow(_Q.row(r)(0),2),
                            std::pow(_Q.row(r)(1),2),
                            std::pow(_Q.row(r)(2),2),
                            _Q.row(r)(0)*_Q.row(r)(1),
                            _Q.row(r)(1)*_Q.row(r)(2),
                            _Q.row(r)(0)*_Q.row(r)(2);
                _Qdelta.row(r) = _qdelta;
            }

            Eigen::MatrixXd _Apq = mass * _P.transpose() * _Qdelta; //3x9
            Eigen::MatrixXd _Aqqinv = (mass * _Qdelta.transpose() * _Qdelta); 
            Eigen::MatrixXd _Aqq = _Aqqinv.inverse();
            //Eigen::MatrixXd _Aqq = pseudoinverse(_Aqqinv, std::numeric_limits<double>::epsilon());//9x9

            Eigen::MatrixXd _A = _Apq * _Aqq; //how do we presearve volume
            Eigen::MatrixXd _R;
            _R.setZero(3,9);
            _R.block(0,0,3,3) = R;
            double beta = 0.7;
            Eigen::MatrixXd T;
            T = beta * _A + (1-beta) * _R; //3x9
            transformedP = (T * _Qdelta.transpose()).transpose();
        }
        Eigen::MatrixXd gt = transformedP.rowwise() + center_of_masst.transpose();
        //flatten gt
        Eigen::VectorXd gt_flatten;
        gt_flatten.resize(gt.rows()*gt.cols());
        Eigen::MatrixXd gtt = gt.transpose();
        gt_flatten = Eigen::Map<Eigen::VectorXd>(gtt.data(), gtt.rows()*gtt.cols());
        
        //update
        // float update_weight = 1.0;
        // if(clusters.size() > 1){
        //     if(ic==0){
        //         update_weight=1.f/2.f;
        //     }else{
        //         update_weight=1.f/(2.f*(clusters.size()-1.f));
        //     }
        // }
        // update_weight=1.f/clusters.size();

        double alpha = 1.0;
        Eigen::VectorXd qdoti_updates = alpha * ((gt_flatten - qi_tmp))/dt;     
        //add to the global qdot at cluster vertices
        for(int iv=0;iv<clusters.at(ic).size();++iv){
            int global_v_idx = clusters.at(ic).at(iv);
            int v_clusters = vertex_cluster_count(global_v_idx);
            float update_weight = 1.f / float(v_clusters);
            //std::cout<<update_weight<<std::endl;
            qdot_holder.segment<3>(global_v_idx*3) += update_weight * (qdot_tmp.segment<3>(global_v_idx*3) + qdoti_updates.segment<3>(iv*3));
            //qdot.segment<3>(clusters.at(ic).at(iv)*3) = update_weight * (qdot_tmp.segment<3>(clusters.at(ic).at(iv)*3) + qdoti_updates.segment<3>(iv*3));
        }
        
        // if(tmp_force.norm() > 60){
        //     std::cout<<"debug"<<std::endl;
        //     std::cout<<"update weight:"<<update_weight<<std::endl;
        //     std::cout<<Eigen::Map<Eigen::MatrixXd>(tmp_force.data(),3,tmp_force.rows()/3).transpose()<<std::endl;
        //     std::cout<<"old:"<<std::endl;
        //     std::cout<<Eigen::Map<Eigen::MatrixXd>(q.data(),3,q.rows()/3).transpose()<<std::endl;
        //     std::cout<<"tmp:"<<std::endl;
        //     std::cout<<Vt.transpose()<<std::endl;
        //     std::cout<<"goals:"<<std::endl;
        //     std::cout<<gt<<std::endl;
        //     std::cout<<"qdot tmp"<<std::endl;
        //     std::cout<<Eigen::Map<Eigen::MatrixXd>(qdot_tmp.data(),3,qdot_tmp.rows()/3).transpose()<<std::endl;
        //     std::cout<<"qdot udpate"<<std::endl;
        //     std::cout<<Eigen::Map<Eigen::MatrixXd>(qdoti_updates.data(),3,qdoti_updates.rows()/3).transpose()<<std::endl;
        //     std::cout<<"qdot"<<std::endl;
        //     std::cout<<Eigen::Map<Eigen::MatrixXd>(qdot.data(),3,qdot.rows()/3).transpose()<<std::endl;
        //     std::cout<<"qdot holder:"<<std::endl;
        //     std::cout<<Eigen::Map<Eigen::MatrixXd>(qdot_holder.data(),3,qdot_holder.rows()/3).transpose()<<std::endl;

        // }
    }
    //update q with qdot
    qdot = qdot_holder;
    if(method==2){
        dt /= 1.2;
    }
    qdot = P.transpose() * P * qdot;
    q = q + dt * qdot;
    //q = P.transpose() * P * q + x0;
}
