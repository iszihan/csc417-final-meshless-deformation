#include <Eigen/Eigen/src/Core/Matrix.h>

class Spatial_hash_fn
{
	int p1 = 73856093;
	int p2 = 19349663;
	int p3 = 83492791;
	int hashtable_size = 0;
	double cellsize = 0.1;

public:
	size_t operator()(const Eigen::Vector3d inputval);
	Spatial_hash_fn(double chosen_cellsize, int chosen_hashtable_size);
};

