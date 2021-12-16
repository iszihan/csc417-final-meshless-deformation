#include "Spatial_hash_fn.h"
#include <math.h>
size_t Spatial_hash_fn::operator()(const Eigen::Vector3d inputval)
{
	int v1, v2, v3, h;
	v1 = (int)floor(inputval(0) / Spatial_hash_fn::cellsize);
	v2 = (int)floor(inputval(1) / Spatial_hash_fn::cellsize);
	v3 = (int)floor(inputval(2) / Spatial_hash_fn::cellsize);
	h = ((v1 * p1) & (v2 * p2) & (v3 * p3)) % hashtable_size;
	return h;

}

Spatial_hash_fn::Spatial_hash_fn(double chosen_cellsize, int chosen_hashtable_size)
{
	cellsize = chosen_cellsize;
	hashtable_size = chosen_hashtable_size;
}
