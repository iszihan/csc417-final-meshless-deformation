#include <collision_detection.h>
#include <iostream>

size_t Spatial_hash_fn::operator()(const Eigen::Vector3d inputval) const
{
    int v1, v2, v3, h;
    v1 = (int)floor(inputval(0) / Spatial_hash_fn::cellsize);
    v2 = (int)floor(inputval(1) / Spatial_hash_fn::cellsize);
    v3 = (int)floor(inputval(2) / Spatial_hash_fn::cellsize);
    h = ((v1 * p1) & (v2 * p2) & (v3 * p3)) % hashtable_size;
    return h;
}
void Spatial_hash_fn::get_mapping(const Eigen::Vector3d inputval, Eigen::Vector3d rounded_values) const
{
    double v1, v2, v3, h;
    rounded_values(0) = floor(inputval(0) / Spatial_hash_fn::cellsize) * Spatial_hash_fn::cellsize;
    rounded_values(1) = floor(inputval(1) / Spatial_hash_fn::cellsize) * Spatial_hash_fn::cellsize;
    rounded_values(2) = floor(inputval(2) / Spatial_hash_fn::cellsize) * Spatial_hash_fn::cellsize;
}


bool point_triangle_intersection(Eigen::Vector3d point, Eigen::VectorXd vertex_list, 
                                 Eigen::RowVector3i face, Eigen::Vector3d & collision_pt)
{
    Eigen::Vector3d v1, v2, v3, phi;
    Eigen::Matrix3d A; 
    v1 = vertex_list.segment<3>(face(0) * 3);
    v2 = vertex_list.segment<3>(face(1) * 3);
    v3 = vertex_list.segment<3>(face(2) * 3);
    Eigen::Vector3d n = ((v2 - v1).cross(v3 - v1)).normalized();
    double dist = (point - v1).dot(n);
    if (abs(dist) > 0.1) {
        return false;
    }
    A.col(0) = v1 - v3;
    A.col(1) = v2 - v3;
    A.col(2) = n;
    Eigen::VectorXd point_projected = point - n * std::max(dist, 0.05);
    phi = A.inverse() * point_projected;
	if (phi(1) >= -0.001 && phi(0) > -0.001 && phi(1) < 1.001 && phi(0) < 1.001 && (phi(1) + phi(0) <= 1.001))
	{
        collision_pt = point_projected;
        return true;
	}
    return false;
}
bool ray_triangle_intersect(
    Eigen::Vector3d ray_origin, Eigen::Vector3d ray_direction, double& t, Eigen::VectorXd vertex_list, Eigen::RowVector3i face)
{
    double min_t = 0.0001;
	
    double local_t = 0;

    
	
    Eigen::Vector3d s = vertex_list.segment<3>(face(0) * 3) - ray_origin;
    Eigen::Vector3d normal;
    Eigen::Matrix3d C;
    C.col(0) = vertex_list.segment<3>(face(0) * 3) - vertex_list.segment<3>(face(1) * 3);
    C.col(1) = vertex_list.segment<3>(face(0) * 3) - vertex_list.segment<3>(face(2) * 3);
    C.col(2) = ray_direction;
    // pre calculating the variable before the computation. The equation from the text book is used.
    double ei_hf, gf_di, dh_eg, ak_jb, jc_al, bl_kc;
    double a = C(0, 0);
    double b = C(1, 0);
    double c = C(2, 0);
    double d = C(0, 1);
    double e = C(1, 1);
    double f = C(2, 1);
    double g = C(0, 2);
    double h = C(1, 2);
    double i = C(2, 2);
    double j = s(0);
    double k = s(1);
    double l = s(2);
    double beta, gamma;
    ei_hf = e * i - h * f;
    gf_di = g * f - d * i;
    dh_eg = d * h - e * g;
    double m = ei_hf * a + gf_di * b + dh_eg * c;
    if (m == 0) {
        return false;
    }
    ak_jb = a * k - j * b;
    jc_al = j * c - a * l;
    bl_kc = b * l - k * c;
    // solve for t, gamma and beta sequentially so the algorithm can terminated earlier if needed
    local_t = -(f * ak_jb + e * jc_al + d * bl_kc) / m;
    if (local_t < min_t) {
        return false;
    }
    gamma = (i * ak_jb + h * jc_al + g * bl_kc) / m;
    if (gamma < 0 || gamma > 1) {
        return false;
    }
    beta = (j * ei_hf + k * gf_di + l * dh_eg) / m;
    if (beta < 0 || beta > 1 - gamma) {
        return false;
    }
    // overwrite t with local value only if the ray actually intersect with the shape
    t = local_t;
	return true;
    ////////////////////////////////////////////////////////////////////////////
}
bool Bounding_box::ray_box_intersection(Eigen::Vector3d origin, Eigen::Vector3d ray_dir)
{
    double txmin, tymin, tzmin, txmax, tymax, tzmax;
    Eigen::Vector3d dir = ray_dir.normalized();
    if (dir(0) > 0) {
        txmax = (max_corner(0, 0) - origin(0)) / dir(0);
        txmin = (min_corner(0, 0) - origin(0)) / dir(0);
    }
    else {
        txmin = (max_corner(0, 0) - origin(0)) / dir(0);
        txmax = (min_corner(0, 0) - origin(0)) / dir(0);
    }
    if (dir(1) > 0) {
        tymax = (max_corner(0, 1) - origin(1)) / dir(1);
        tymin = (min_corner(0, 1) - origin(1)) / dir(1);
    }
    else {
        tymin = (max_corner(0, 1) - origin(1)) / dir(1);
        tymax = (min_corner(0, 1) - origin(1)) / dir(1);
    }
    if (dir(2) > 0) {
        tzmax = (max_corner(0, 2) - origin(2)) / dir(2);
        tzmin = (min_corner(0, 2) - origin(2)) / dir(2);
    }
    else {
        tzmin = (max_corner(0, 2) - origin(2)) / dir(2);
        tzmax = (min_corner(0, 2) - origin(2)) / dir(2);
    }
    if (std::max(std::max(txmin, tymin), tzmin) < std::min(std::min(txmax, tymax), tzmax)) {
        return true;
    }
    else if (std::max(std::max(txmin, tymin), tzmin) < 0 && std::min(std::min(txmax, tymax), tzmax) > 0) {
        return true;
    }
    return false;
}
bool Bounding_box::add_face(int face_id)
{
    face_list.push_back(face_id);
    return true;
}
bool Bounding_box::get_face_list(std::vector<int>& rtv_face_list)
{
	if (face_list.size() == 0)
	{
        return false;
	} else
	{
        rtv_face_list = face_list;
        return true;
	}
}
Bounding_box::Bounding_box(Eigen::Vector3d min, Eigen::Vector3d max)
{
    min_corner = min;
    max_corner = max;
}

void collision_detection(std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>> &collisions,
                         unsigned int moving_obj_geometry_idx,
                         unsigned int still_obj_geometry_idx,
                         scene_object obj1, scene_object obj2){       

    Eigen::VectorXd q = std::get<8>(obj1);
    Eigen::VectorXd q2 = std::get<8>(obj2);
    Eigen::MatrixXd sV = std::get<1>(obj2);
    Eigen::MatrixXi sF = std::get<2>(obj2);
    if (std::get<0>(obj2) == 0) {
        for (int vi = 0; vi < q.rows() / 3; ++vi) {
            Eigen::Vector3d curr_v = q.segment<3>(3 * vi);
            //if it is plane 
            Eigen::Vector3d pos = sV.row(0);
            Eigen::Vector3d e10 = sV.row(sF.row(0)(1)) - sV.row(sF.row(0)(0));
            Eigen::Vector3d e20 = sV.row(sF.row(0)(2)) - sV.row(sF.row(0)(0));
            Eigen::Vector3d dir = e10.cross(e20);
            dir = dir.normalized();
            // std::cout<<"plane point"<<std::endl;
            // std::cout<<pos<<std::endl;
            // std::cout<<"plane normal"<<std::endl;
            // std::cout<<dir<<std::endl;

            double dist = (curr_v - pos).dot(dir);
            if (dist < 0.2) {
                std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int> collision_info;
                std::get<0>(collision_info) = curr_v - dir * dist;
                std::get<1>(collision_info) = -dir; // negative sign is added to make it the vertex normal
                std::get<2>(collision_info) = vi;
                std::get<3>(collision_info) = sF.row(0)(1);
                std::get<4>(collision_info) = still_obj_geometry_idx;
                collisions.push_back(collision_info);
            }
        }
    } else if (std::get<0>(obj2) == -1)
    {
        for (int vi = 0; vi < q.rows() / 3; ++vi) {
            Eigen::Vector3d curr_v = q.segment<3>(3 * vi);
            //if it is plane 
            Eigen::Vector3d pos = sV.row(0);
            Eigen::Vector3d e10 = sV.row(sF.row(0)(1)) - sV.row(sF.row(0)(0));
            Eigen::Vector3d e20 = sV.row(sF.row(0)(2)) - sV.row(sF.row(0)(0));
            Eigen::Vector3d dir = e10.cross(e20);
            dir = dir.normalized();

            double dist = (curr_v(1) - pos(1));
            if (dist <= 0.2) {
                std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int> collision_info;
                std::get<0>(collision_info) = curr_v - dir * dist;
                std::get<1>(collision_info) = -dir; // negative sign is added to make it the vertex normal
                std::get<2>(collision_info) = vi;
                std::get<3>(collision_info) = sF.row(0)(1);
                std::get<4>(collision_info) = still_obj_geometry_idx;
                collisions.push_back(collision_info);
            }
        }
    } else
    {
        std::vector<std::vector<int>> VtF = std::get<17>(obj1);
    	Eigen::MatrixXi F1 = std::get<2>(obj1);
        for (int vi = 0; vi < q.rows() / 3; vi++) {
            bool found_p2, found;
            found = false;
            int p2_index;
            found_p2 = false;
            Eigen::Vector3d p2, n1, p3, n2;
            double faces_count = (double)VtF.at(vi).size();
            for (int fi_vi = 0; fi_vi < faces_count; fi_vi++)
            {
                Eigen::Vector3d v0, v1, v2, n_of_vi;
            	v0 = q.segment<3>(F1.row(VtF.at(vi).at(fi_vi))(0) * 3);
                v1 = q.segment<3>(F1.row(VtF.at(vi).at(fi_vi))(1) * 3);
                v2 = q.segment<3>(F1.row(VtF.at(vi).at(fi_vi))(2) * 3);
                n_of_vi = (v1 - v0).cross(v2 - v0);
                n1 = n1 + 1.0 / faces_count * n_of_vi;
            }
            n1 = n1.normalized();
        	// see if there are any intersections between the vertex p1 on obj1 and faces on obj2
            for (int fi = 0; fi < sF.rows(); fi++)
            {
                found = false;
                Eigen::Vector3d intersection_pt;
                double t;
            	// obtain the vecters normal (using an average of the face normals around the vertex)
                if (found_p2) {
                    found = ray_triangle_intersect(p2, n1, t, q2, sF.row(fi));
                } else
                {
                    found = ray_triangle_intersect(q.segment<3>(vi * 3), -n1, t, q2, sF.row(fi));
                }
            	// if there is on point of intersection, 
            	if (found && !found_p2)
            	{
                    p2 = q.segment<3>(vi * 3) + t * (-n1);
                    Eigen::Vector3d v0, v1, v2;
                    v0 = q2.segment<3>(sF.row(fi)(0) * 3);
                    v1 = q2.segment<3>(sF.row(fi)(1) * 3);
                    v2 = q2.segment<3>(sF.row(fi)(2) * 3);
                    n2 = (v1 - v0).cross(v2 - v0);
                    n2 = n2.normalized();
                    if (n1.dot(n2) > 0)
                    {
                        break;
                    }
                    found_p2 = true;
                    p2_index = sF.row(fi)(0);
            	} else if (found && found_p2)
            	{
            		p3 = q.segment<3>(vi * 3) + t * (n1);
            		if ((p2 - q.segment<3>(vi * 3)).norm() <= (p3 - q.segment<3>(vi * 3)).norm())
            		{
                        std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int> collision_info;
                        std::get<0>(collision_info) = p2;
                        std::get<1>(collision_info) = -n2;
                        std::get<2>(collision_info) = vi;
                        std::get<3>(collision_info) = p2_index;
                        std::get<4>(collision_info) = still_obj_geometry_idx;
                        collisions.push_back(collision_info);
                        break;
            		}
            	}
            }
        }
    }
}
void compute_vertex_face_list(Eigen::MatrixXd V, Eigen::MatrixXi F, std::vector<std::vector<int>> &V2F)
{
    V2F.clear();
	for (int i = 0; i < V.rows(); i++)
	{
        std::vector<int> current_vertex;
        current_vertex.clear();
		for (int j = 0; j < F.rows(); j++)
		{
            for (int v = 0; v < 3; v++)
            {
	            if (F.row(j)(v) == i)
	            {
                    current_vertex.push_back(j);
                    break;
	            }
            }
		}
        V2F.push_back(current_vertex);
	}
}

bool precomputation(scene_object obj1, scene_object obj2)
{
    Eigen::VectorXd q1 = std::get<8>(obj1);
    Eigen::VectorXd q2 = std::get<8>(obj2);
    Eigen::Vector3d com1, com2;
    double rad1, rad2;
	
	
    if (std::get<0>(obj2) >= 1){
        com1 = std::get<16>(obj1);
        rad1 = std::get<15>(obj1);
        com2 = std::get<16>(obj2);
        rad2 = std::get<15>(obj2);
        return ((com2 - com1).norm() + 0.1 <= std::min(rad1, rad2));
    } else
    {
        com1 = std::get<16>(obj1);
        rad1 = std::get<15>(obj1);
    	Eigen::MatrixXd sV = std::get<1>(obj2);
        Eigen::MatrixXi sF = std::get<2>(obj2);

        Eigen::Vector3d pos = sV.row(sF.row(0)(0)).transpose();
        Eigen::Vector3d e10 = sV.row(sF.row(0)(1)).transpose() - sV.row(sF.row(0)(0)).transpose();
        Eigen::Vector3d e20 = sV.row(sF.row(0)(2)).transpose() - sV.row(sF.row(0)(0)).transpose();
        Eigen::Vector3d dir = e10.cross(e20);
        dir = dir.normalized();
        
        double dist = abs((com1 - pos).dot(dir));
    	return dist - rad1 <= 0.1;
    }
}

void construct_spatial_hash_table(scene_object & obj)
{
    std::unordered_map <Eigen::Vector3d, std::vector<int>, Spatial_hash_fn> ht;
    std::vector<int> occupied_hash_keys;
    Spatial_hash_fn hf;
    occupied_hash_keys.clear();
    ht.clear();
    Eigen::MatrixXi F = std::get<2>(obj);
    Eigen::VectorXd q = std::get<8>(obj);
    Eigen::Vector3d v0, v1, v2, max_cor, min_cor;
    for (int i = 0; i < F.rows(); i++)
    {
        v0 = q.segment<3>(F.row(i)(0) * 3);
        v1 = q.segment<3>(F.row(i)(1) * 3);
        v2 = q.segment<3>(F.row(i)(2) * 3);
        max_cor(0) = std::max(std::max(v0(0), v1(0)), v2(0));
        max_cor(1) = std::max(std::max(v0(1), v1(1)), v2(1));
        max_cor(2) = std::max(std::max(v0(2), v1(2)), v2(2));
        min_cor(0) = std::min(std::min(v0(0), v1(0)), v2(0));
        min_cor(1) = std::min(std::min(v0(1), v1(1)), v2(1));
        min_cor(2) = std::min(std::min(v0(2), v1(2)), v2(2));
    }
    return;
}
