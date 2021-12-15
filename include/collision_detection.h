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

//typedef std::tuple<int, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::SparseMatrixd,
//    Eigen::SparseMatrixd, Eigen::Vector3d, Eigen::VectorXd, Eigen::VectorXd, Eigen::SparseMatrixd, Eigen::VectorXd, double, Eigen::Vector3d> scene_object;

typedef std::tuple<int,                           //moving or still
    Eigen::MatrixXd,               //V
    Eigen::MatrixXi,               //F
    Eigen::MatrixXd,               //V_skin
    Eigen::MatrixXi,               //F_skin
    Eigen::SparseMatrixd,          //N skinning matrix
    Eigen::SparseMatrixd,          //M
    Eigen::Vector3d,               //center of mass
    Eigen::VectorXd,               //q
    Eigen::VectorXd,               //qdot
    Eigen::SparseMatrixd,          //P
    Eigen::VectorXd,               //x0
    Eigen::VectorXd,               //gravity
    std::vector<std::vector<int>>, //clusters
    Eigen::MatrixXd,               //Q = V-center_of_mass
    double,                         //distance to com
    Eigen::Vector3d,                 // com that moves with time
    std::vector<std::vector<int>>    // vertex face list
> scene_object;

void collision_detection(std::vector<std::pair<Eigen::Vector3d, unsigned int>> &collisions,
                         unsigned int moving_obj_type_id,
                         unsigned int still_obj_type_id,
                         scene_object obj1, scene_object obj2);
bool precomputation(scene_object obj1, scene_object obj2);
        