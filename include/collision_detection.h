#include <Eigen/Dense>
#include <EigenTypes.h>
#include <vector>
#include <tuple>
typedef std::tuple<int,                           //0 moving or still
                   Eigen::MatrixXd,               //1 V
                   Eigen::MatrixXi,               //2 F
                   Eigen::MatrixXd,               //3 V_skin
                   Eigen::MatrixXi,               //4 F_skin
                   Eigen::SparseMatrixd,          //5 N skinning matrix
                   Eigen::SparseMatrixd,          //6 M -- this is not actually used anywhere?
                   std::vector<Eigen::Vector3d>,  //7 center of mass for clusters
                   Eigen::VectorXd,               //8 q
                   Eigen::VectorXd,               //9 qdot
                   Eigen::SparseMatrixd,          //10 P
                   Eigen::VectorXd,               //11 x0
                   Eigen::VectorXd,               //12 gravity
                   std::vector<std::vector<int>>, //13 clusters of pair of vertex indices and positions
                   std::vector<Eigen::MatrixXd>,  //14 Q = V-center_of_mass for clusters
                   double,                        //15 distance to com
                   Eigen::Vector3d                //16 com that moves with time
                   >scene_object;
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
                         unsigned int moving_obj_type_id,
                         unsigned int still_obj_type_id,
                         scene_object obj1, scene_object obj2);

bool precomputation(scene_object obj1, scene_object obj2);
        