#include <collision_detection.h>
#include <iostream>

void collision_detection(std::vector<std::pair<Eigen::Vector3d, unsigned int>> &collisions,
                         unsigned int moving_obj_id,
                         unsigned int still_obj_id,
                         Eigen::Ref<Eigen::VectorXd> q, 
                         Eigen::Ref<Eigen::MatrixXd> sV, 
                         Eigen::Ref<Eigen::MatrixXi> sF){
    //TODO: need to handle different object differently? 
    //PLANE: check against the sides of the plane, just need plane normal and a plane vertex 
    //OTHERS: check against every vertices and compute distance -- no inside/outside check here? not needed?
    for(int vi=0;vi<q.rows()/3;++vi){
        Eigen::Vector3d curr_v = q.segment<3>(3*vi);
        //if it is plane 
        Eigen::Vector3d pos = sV.row(0);
        Eigen::Vector3d e10 = sV.row(sF.row(0)(1))-sV.row(sF.row(0)(0));
        Eigen::Vector3d e20 = sV.row(sF.row(0)(2))-sV.row(sF.row(0)(0));
        Eigen::Vector3d dir = e10.cross(e20);
        dir = dir.normalized();
        // std::cout<<"plane point"<<std::endl;
        // std::cout<<pos<<std::endl;
        // std::cout<<"plane normal"<<std::endl;
        // std::cout<<dir<<std::endl;
        
        double dist = (curr_v - pos).dot(dir);
        if(dist < 0.2){
            //save collision
            //std::cout<<"dist:"<<dist<<std::endl;
            //std::cout<<"current v"<<std::endl;
            //std::cout<<curr_v<<std::endl;
            //std::cout<<"plane v"<<std::endl;
            //std::cout<<curr_v - dir * dist<<std::endl;
            collisions.push_back(std::make_pair(curr_v-dir*dist, vi));
        }
    }
}

bool precomputation(scene_object obj1, scene_object obj2)
{
    Eigen::VectorXd q1 = std::get<8>(obj1);
    Eigen::VectorXd q2 = std::get<8>(obj2);
    Eigen::Vector3d com1, com2;
    double rad1, rad2;
	
	
    if (std::get<0>(obj2) >= 1){
        com1 = std::get<13>(obj1);
        rad1 = std::get<12>(obj1);
        com2 = std::get<13>(obj2);
        rad2 = std::get<12>(obj2);
        return ((com2 - com1).norm() + 0.1 <= std::min(rad1, rad2));
    } else
    {
        com1 = std::get<13>(obj1);
        rad1 = std::get<12>(obj1);
    	Eigen::MatrixXd sV = std::get<1>(obj2);
        Eigen::MatrixXi sF = std::get<2>(obj2);

        Eigen::Vector3d pos = sV.row(0).transpose();
        Eigen::Vector3d e10 = sV.row(sF.row(0)(1)).transpose() - sV.row(sF.row(0)(0)).transpose();
        Eigen::Vector3d e20 = sV.row(sF.row(0)(2)).transpose() - sV.row(sF.row(0)(0)).transpose();
        Eigen::Vector3d dir = e10.cross(e20);
        dir = dir.normalized();
        double dist = abs((com1 - pos).dot(dir));
        return dist - rad1 <= 0.1;
    }
}
