#include <collision_box_floor.h>
#include <iostream>

//Input:
//  R - rotation matrix for rigid body
//  p - world space position of center-of-mass
//  obj_id - id of object being checked for collision with the floor
//  V - the nx3 matrix of undeformed vertex positions [n x 3]
//  dir - the outward facing normal for the floor
//  pos - the world space position of a point on the floor plane
//Output:
//  x - world space, per-vertex collision points. Computed as any vertex that is on the "wrong side" of the floor plane
//  n - collision normals, one for each collision point. These point away from the floor. 
//  objs - Pairs of ids for objects involved in collisions. The first id is for the object, away from which the normal points. The second id 
//  is for the object towards which the normal points. The floor has an id of -1. 
void collision_box_floor(std::vector<Eigen::Vector3d> &x, std::vector<Eigen::Vector3d> &n, std::vector<std::pair<int,int>> &objs,
                         Eigen::Ref<const Eigen::Matrix3d> R, Eigen::Ref<const Eigen::Vector3d> p, unsigned int obj_id,
                         Eigen::Ref<Eigen::MatrixXd> V, 
                         Eigen::Ref<const Eigen::Vector3d> dir, Eigen::Ref<const Eigen::Vector3d> pos) {

    // get deformed vertex positions
    Eigen::MatrixXd ptiled = p.replicate(1,V.rows());
    Eigen::MatrixXd v = (R * V.transpose() + ptiled).transpose();
    // for each vertex, check whether it is on the wrong side of the floor plane 
    for(int i=0;i<v.rows();++i){
        Eigen::Vector3d curr_v = v.row(i);
        double sd2floor = (curr_v - pos).dot(dir);
        if(sd2floor<0){
            //save collision info
            x.push_back(curr_v);
            n.push_back(dir);
            objs.push_back(std::make_pair(-1, obj_id));
        }
    }
}