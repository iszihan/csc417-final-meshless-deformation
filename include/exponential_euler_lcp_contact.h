#include <Eigen/Dense>
#include <EigenTypes.h>
#include <rodrigues.h>
#include <iostream>
#include <rigid_body_jacobian.h>
#include <inverse_rigid_body.h>

//Input:
//  q - 12n vector where n is the number of rigid bodies. Each rigid body is stored as 12 doubles. 
//      The first 9 doubles are the columns of a 3x3 rotation matrix and the final 3 doubles are the world space position of the object's center of mass.
//  qdot - 6n vector of generalied velocities. The first 3 doubles of each body are the world space angular velocity and 
//         the second 3 are the world space linear velocity.
//  dt - the integration time step
//  masses - a vector to mass matrices for each rigid body
//  forces - a 6n vector of generalized forces for n rigid bodies. The first 3 doubles of each rigid body are the torques acting on the object
//           while the second 3 doubles are the linear forces.
//  n - list of collision normals
//  x - list of world space collision points
//  obj - list of collision object ids 
//Output:
//  q - updated generalized coordinates 
//  qdot - updated generalized velocities 
inline void exponential_euler_lcp_contact(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                                          std::vector<Eigen::Matrix66d> &masses, Eigen::Ref<const Eigen::VectorXd> forces,
                                          std::vector<Eigen::Vector3d> &n, std::vector<Eigen::Vector3d> &x, std::vector<std::pair<int,int> > &obj) {

    int ncontacts = n.size();
    int maxitr = 10;
    int it = 0;

    //compute qdot_star without constraints 
    // Eigen::Vector6d qdot_star;
    // Eigen::Matrix3d inertia = masses[0].block(0,0,3,3);
    // Eigen::Matrix3d Rt;
    // Rt<<q(0),q(3),q(6),
    //     q(1),q(4),q(7),
    //     q(2),q(5),q(8);
    // Eigen::Matrix3d RIRT = Rt * inertia * Rt.transpose();
    // Eigen::Matrix3d exp_wt;
    // Eigen::Vector3d omega = qdot.head(3);
    // rodrigues(exp_wt,omega*dt);
    // Eigen::Vector3d b = RIRT * omega + dt * omega.cross(RIRT * omega) + dt * forces.head(3);
    // qdot_star.head(3) = RIRT.inverse() * b; 
    // qdot_star.tail(3) = (masses[0](3,3)*qdot.tail(3)+dt*forces.tail(3)) / masses[0](3,3);

    //compute qdot_star with M matrix constructed
    Eigen::Matrix66d M;
    M.setZero(6,6);
    Eigen::Matrix3d inertia = masses[0].block(0,0,3,3);
    Eigen::Matrix3d Rt;
    Rt<<q(0),q(3),q(6),
        q(1),q(4),q(7),
        q(2),q(5),q(8);    
    Eigen::Matrix3d RIRT = Rt * inertia * Rt.transpose();
    M.block(0,0,3,3) = RIRT;
    M.block(3,3,3,3) = masses[0](3,3) * Eigen::Matrix3d::Identity();
    
    Eigen::Vector6d f;
    Eigen::Matrix3d exp_wt;
    Eigen::Vector3d omega = qdot.head(3);
    rodrigues(exp_wt,omega*dt);
    f.head(3) = (omega.cross(RIRT*omega) + forces.head(3));
    f.tail(3) = forces.tail(3);

    Eigen::Vector6d qdot_star = qdot + dt * M.inverse() * f;

    //projected gauss-seidel to find alphas
    Eigen::VectorXd alpha(ncontacts,1);
    alpha.setZero(ncontacts,1);
    while(it < maxitr){
        for(int ic=0; ic<ncontacts; ic++){
            Eigen::Vector3d n_contact = n[ic];
            Eigen::Vector3d x_contact = x[ic];
            Eigen::Vector3d X_contact;
            inverse_rigid_body(X_contact, x_contact, Rt, q.tail(3));
            Eigen::Matrix36d J;
            rigid_body_jacobian(J, Rt, q.tail(3), X_contact);
            // delta term deals with contact point only
            double delta_ic = dt * (n_contact.transpose() * J * M.inverse() * J.transpose() * n_contact)(0);
            // gamma terms deal with all the other contact points
            Eigen::Vector6d b;
            b.setZero(6,1);
            for(int i= 0; i<ncontacts; i++){
                if(i!=ic){
                    Eigen::Vector3d ni_contact = n[i];
                    Eigen::Vector3d xi_contact = x[i];
                    Eigen::Vector3d Xi_contact;
                    inverse_rigid_body(Xi_contact, xi_contact, Rt, q.tail(3));
                    Eigen::Matrix36d Ji;
                    rigid_body_jacobian(Ji, Rt, q.tail(3), Xi_contact);
                    b += dt * M.inverse() * Ji.transpose() * ni_contact * alpha(i);
                }
            }
            double gamma_ic = n_contact.transpose() * J * (qdot_star + b);
            alpha(ic) = std::max(0.0,-delta_ic/gamma_ic);
        }
        it++;
    }

    //compute qdott_1 with the alphas
    for(int i= 0; i<ncontacts; i++){
        Eigen::Vector3d ni_contact = n[i];
        Eigen::Vector3d xi_contact = x[i];
        Eigen::Vector3d Xi_contact;
        inverse_rigid_body(Xi_contact, xi_contact, Rt, q.tail(3));
        Eigen::Matrix36d Ji;
        rigid_body_jacobian(Ji, Rt, q.tail(3), Xi_contact);
        qdot_star += dt * M.inverse() * Ji.transpose() * ni_contact * alpha(i);
    }

    qdot = qdot_star;
    q.tail(3) = q.tail(3) + dt * qdot.tail(3);
    Eigen::Matrix3d Rt1 = exp_wt * Rt;
    q(0) = Rt1(0,0);
    q(1) = Rt1(1,0);
    q(2) = Rt1(2,0);
    q(3) = Rt1(0,1);
    q(4) = Rt1(1,1);
    q(5) = Rt1(2,1);
    q(6) = Rt1(0,2);
    q(7) = Rt1(1,2);
    q(8) = Rt1(2,2);
}