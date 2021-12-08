#include <inertia_matrix.h>
#include <cassert>

//compute inertia matrix and volume by integrating on surfaces
void inertia_matrix(Eigen::Matrix3d &I, Eigen::Vector3d & center, double &mass, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double density) {
    mass = 0.0;
    I.setZero(3,3);
    center.setZero();
    for(int i=0; i<F.rows();++i){
        Eigen::Vector3d X0 = V.row(F.row(i)(0));
        Eigen::Vector3d X1 = V.row(F.row(i)(1));
        Eigen::Vector3d X2 = V.row(F.row(i)(2));
        Eigen::Vector3d crp = (X1-X0).cross(X2-X0);
        Eigen::Vector3d normal = crp.normalized();
        double area = crp.norm() / 2.0; 

        //compute mass
        double curr_mass = density * area * (X0(0) + X1(0) + X2(0)) * normal(0) / 3.0;
        mass += curr_mass;

        Eigen::Vector3d int_com;
        int_com(0) = normal(0) * ((X0(0)*X1(0))/1.2E+1+(X0(0)*X2(0))/1.2E+1+(X1(0)*X2(0))/1.2E+1+(X0(0)*X0(0))/1.2E+1+(X1(0)*X1(0))/1.2E+1+(X2(0)*X2(0))/1.2E+1);
        int_com(1) = normal(1) * ((X0(1)*X1(1))/1.2E+1+(X0(1)*X2(1))/1.2E+1+(X1(1)*X2(1))/1.2E+1+(X0(1)*X0(1))/1.2E+1+(X1(1)*X1(1))/1.2E+1+(X2(1)*X2(1))/1.2E+1);
        int_com(2) = normal(2) * ((X0(2)*X1(2))/1.2E+1+(X0(2)*X2(2))/1.2E+1+(X1(2)*X2(2))/1.2E+1+(X0(2)*X0(2))/1.2E+1+(X1(2)*X1(2))/1.2E+1+(X2(2)*X2(2))/1.2E+1);
        Eigen::Vector3d curr_com = density * area * int_com;
        center += curr_com;
    }
    center /= mass;
    //std::cout<<"center of mass:"<<center<<"| mass:"<<mass<<std::endl;

    for(int i=0; i<F.rows();++i){
        Eigen::Vector3d X0 = V.row(F.row(i)(0)).transpose()- center;
        Eigen::Vector3d X1 = V.row(F.row(i)(1)).transpose() - center;
        Eigen::Vector3d X2 = V.row(F.row(i)(2)).transpose() - center;
        Eigen::Vector3d crp = (X1-X0).cross(X2-X0);
        Eigen::Vector3d normal = crp.normalized();
        double area = crp.norm() / 2.0; 

        //compute inertia 
        Eigen::Matrix3d curr_I;
        double int_xx2 = (X0(0)*(X1(0)*X1(0)))/2.0E+1+((X0(0)*X0(0))*X1(0))/2.0E+1+(X0(0)*(X2(0)*X2(0)))/2.0E+1+((X0(0)*X0(0))*X2(0))/2.0E+1+(X1(0)*(X2(0)*X2(0)))/2.0E+1+((X1(0)*X1(0))*X2(0))/2.0E+1+(X0(0)*X0(0)*X0(0))/2.0E+1+(X1(0)*X1(0)*X1(0))/2.0E+1+(X2(0)*X2(0)*X2(0))/2.0E+1+(X0(0)*X1(0)*X2(0))/2.0E+1;
        double int_xy2 = (X0(1)*(X1(1)*X1(1)))/2.0E+1+((X0(1)*X0(1))*X1(1))/2.0E+1+(X0(1)*(X2(1)*X2(1)))/2.0E+1+((X0(1)*X0(1))*X2(1))/2.0E+1+(X1(1)*(X2(1)*X2(1)))/2.0E+1+((X1(1)*X1(1))*X2(1))/2.0E+1+(X0(1)*X0(1)*X0(1))/2.0E+1+(X1(1)*X1(1)*X1(1))/2.0E+1+(X2(1)*X2(1)*X2(1))/2.0E+1+(X0(1)*X1(1)*X2(1))/2.0E+1;
        double int_xz2 = (X0(2)*(X1(2)*X1(2)))/2.0E+1+((X0(2)*X0(2))*X1(2))/2.0E+1+(X0(2)*(X2(2)*X2(2)))/2.0E+1+((X0(2)*X0(2))*X2(2))/2.0E+1+(X1(2)*(X2(2)*X2(2)))/2.0E+1+((X1(2)*X1(2))*X2(2))/2.0E+1+(X0(2)*X0(2)*X0(2))/2.0E+1+(X1(2)*X1(2)*X1(2))/2.0E+1+(X2(2)*X2(2)*X2(2))/2.0E+1+(X0(2)*X1(2)*X2(2))/2.0E+1;
        double int_xxxy = ((X0(0)*X0(0))*X0(1))/2.0E+1+((X0(0)*X0(0))*X1(1))/6.0E+1+((X1(0)*X1(0))*X0(1))/6.0E+1+((X0(0)*X0(0))*X2(1))/6.0E+1+((X1(0)*X1(0))*X1(1))/2.0E+1+((X2(0)*X2(0))*X0(1))/6.0E+1+((X1(0)*X1(0))*X2(1))/6.0E+1+((X2(0)*X2(0))*X1(1))/6.0E+1+((X2(0)*X2(0))*X2(1))/2.0E+1+(X0(0)*X1(0)*X0(1))/3.0E+1+(X0(0)*X1(0)*X1(1))/3.0E+1+(X0(0)*X2(0)*X0(1))/3.0E+1+(X0(0)*X1(0)*X2(1))/6.0E+1+(X0(0)*X2(0)*X1(1))/6.0E+1+(X1(0)*X2(0)*X0(1))/6.0E+1+(X0(0)*X2(0)*X2(1))/3.0E+1+(X1(0)*X2(0)*X1(1))/3.0E+1+(X1(0)*X2(0)*X2(1))/3.0E+1;
        double int_xxxz = ((X0(0)*X0(0))*X0(2))/2.0E+1+((X0(0)*X0(0))*X1(2))/6.0E+1+((X1(0)*X1(0))*X0(2))/6.0E+1+((X0(0)*X0(0))*X2(2))/6.0E+1+((X1(0)*X1(0))*X1(2))/2.0E+1+((X2(0)*X2(0))*X0(2))/6.0E+1+((X1(0)*X1(0))*X2(2))/6.0E+1+((X2(0)*X2(0))*X1(2))/6.0E+1+((X2(0)*X2(0))*X2(2))/2.0E+1+(X0(0)*X1(0)*X0(2))/3.0E+1+(X0(0)*X1(0)*X1(2))/3.0E+1+(X0(0)*X2(0)*X0(2))/3.0E+1+(X0(0)*X1(0)*X2(2))/6.0E+1+(X0(0)*X2(0)*X1(2))/6.0E+1+(X1(0)*X2(0)*X0(2))/6.0E+1+(X0(0)*X2(0)*X2(2))/3.0E+1+(X1(0)*X2(0)*X1(2))/3.0E+1+(X1(0)*X2(0)*X2(2))/3.0E+1;
        double int_xyxz = (X0(1)*(X1(1)*X1(1)))/2.0E+1+((X0(1)*X0(1))*X1(1))/2.0E+1+(X0(1)*(X2(1)*X2(1)))/2.0E+1+((X0(1)*X0(1))*X2(1))/2.0E+1+(X1(1)*(X2(1)*X2(1)))/2.0E+1+((X1(1)*X1(1))*X2(1))/2.0E+1+(X0(1)*X0(1)*X0(1))/2.0E+1+(X1(1)*X1(1)*X1(1))/2.0E+1+(X2(1)*X2(1)*X2(1))/2.0E+1+(X0(1)*X1(1)*X2(1))/2.0E+1;
        curr_I(0,0) = 2.0 * density * area * ( normal(1) * int_xy2 + normal(2) * int_xz2) / 3.0;
        curr_I(1,1) = 2.0 * density * area * ( normal(0) * int_xx2 + normal(2) * int_xz2) / 3.0;
        curr_I(2,2) = 2.0 * density * area * ( normal(0) * int_xx2 + normal(1) * int_xy2) / 3.0;
        curr_I(0,1) = - density * area * normal(0) * int_xxxy;
        curr_I(1,0) = curr_I(0,1);
        curr_I(0,2) = - density * area * normal(0) * int_xxxz;
        curr_I(2,0) = curr_I(0,2);
        curr_I(1,2) = - density * area * normal(0) * int_xyxz;
        curr_I(2,1) = curr_I(1,2);
        I += curr_I;
    }

    
}