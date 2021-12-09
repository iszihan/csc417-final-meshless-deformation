#include <Eigen/Dense>
#include <EigenTypes.h>
#include <vector>
#include <tuple>

//Input:
//  R - rotation matrix for rigid body
//  p - world space position of center-of-mass
//  obj_id - id of object being checked for collision with the floor
//  V - the nx3 matrix of undeformed vertex positions
//  dir - the outward facing normal for the floor
//  pos - the world space position of a point on the floor plane
//Output:
//  x - world space, per-vertex collision points. Computed as any vertex that is on the "wrong side" of the floor plane
//  n - collision normals, one for each collision point. These point away from the floor. 
//  objs - Pairs of ids for objects involved in collisions. The first id is for the object, away from which the normal points. The second id 
//  is for the object towards which the normal points. The floor has an id of -1. 
void collision_detection(std::vector<std::pair<Eigen::Vector3d, unsigned int>> &collisions,
                         unsigned int moving_obj_id,
                         unsigned int still_obj_id,
                         Eigen::Ref<Eigen::VectorXd> mV, 
                         Eigen::Ref<Eigen::MatrixXd> sV, 
                         Eigen::Ref<Eigen::MatrixXi> sf);
        