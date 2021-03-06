#include <Eigen/Dense>
#include <EigenTypes.h>
#include <vector>
#include <unordered_map>
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


class Spatial_hash_fn
{
	int p1 = 73856093;
	int p2 = 19349663;
	int p3 = 83492791;
	int hashtable_size = 500;
public:
	double cellsize = 0.45;

public:
	size_t operator()(const Eigen::Vector3d inputval) const;
	void get_mapping(const Eigen::Vector3d inputval, Eigen::Vector3d & rounded_values) const;
};

class Bounding_box
{
	std::vector<int> face_list;
	Eigen::Vector3d min_corner;
	Eigen::Vector3d max_corner;
public:
	Bounding_box(Eigen::Vector3d min, Eigen::Vector3d max);
	Bounding_box();
	bool ray_box_intersection(Eigen::Vector3d origin, Eigen::Vector3d ray_dir);
	bool add_face(int face_id);
	bool get_face_list(std::vector<int> & rtv_face_list);

};

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
	Eigen::Vector3d,                //16 com that moves with time
	std::vector<std::vector<int>>,   //17  vertex face list
	std::unordered_map <int, Bounding_box>, // 18 spacial hash table
	std::vector<int>,				// 19 occupied list
	Eigen::VectorXi,				// 20 per-vertex cluster counts
	std::vector<Eigen::Matrix3d>	// 21 per-cluster Sp for plasticity
> scene_object;

void collision_detection_old(std::vector<std::pair<Eigen::Vector3d, unsigned int>> &collisions,
                             unsigned int moving_obj_type_id,
                             unsigned int still_obj_type_id,
						     scene_object obj1, scene_object obj2);
void collision_detection(std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>> &collisions,
                         unsigned int moving_obj_type_id,
                         unsigned int still_obj_type_id,
                         scene_object obj1, scene_object obj2);

void collision_detection_with_optimization(std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>>& collisions,
	unsigned int moving_obj_type_id,
	unsigned int still_obj_type_id,
	scene_object obj1, scene_object obj2);

void compute_vertex_face_list(Eigen::MatrixXd V, Eigen::MatrixXi F, std::vector<std::vector<int>>& V2F);

bool precomputation(scene_object obj1, scene_object obj2);

void construct_spatial_hash_table(scene_object & obj);


        