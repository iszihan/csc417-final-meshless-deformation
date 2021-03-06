#include <collision_detection.h>
#include <iostream>

size_t Spatial_hash_fn::operator()(const Eigen::Vector3d inputval) const
{
    int v1, v2, v3;
	unsigned int h;
    v1 = floor(inputval(0) / Spatial_hash_fn::cellsize);
    v2 = floor(inputval(1) / Spatial_hash_fn::cellsize);
    v3 = floor(inputval(2) / Spatial_hash_fn::cellsize);
    h = (unsigned int) ((v1 * p1) & (v2 * p2) & (v3 * p3)) % hashtable_size;
    return h;
}
void Spatial_hash_fn::get_mapping(const Eigen::Vector3d inputval, Eigen::Vector3d & rounded_values) const
{
    double v1, v2, v3;
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
        txmax = (max_corner(0) - origin(0)) / dir(0);
        txmin = (min_corner(0) - origin(0)) / dir(0);
    }
    else {
        txmin = (max_corner(0) - origin(0)) / dir(0);
        txmax = (min_corner(0) - origin(0)) / dir(0);
    }
    if (dir(1) > 0) {
        tymax = (max_corner(1) - origin(1)) / dir(1);
        tymin = (min_corner(1) - origin(1)) / dir(1);
    }
    else {
        tymin = (max_corner(1) - origin(1)) / dir(1);
        tymax = (min_corner(1) - origin(1)) / dir(1);
    }
    if (dir(2) > 0) {
        tzmax = (max_corner(2) - origin(2)) / dir(2);
        tzmin = (min_corner(2) - origin(2)) / dir(2);
    }
    else {
        tzmin = (max_corner(2) - origin(2)) / dir(2);
        tzmax = (min_corner(2) - origin(2)) / dir(2);
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
    Bounding_box::min_corner = min;
    Bounding_box::max_corner = max;
}
Bounding_box::Bounding_box() {}

void collision_detection_old(std::vector<std::pair<Eigen::Vector3d, unsigned int>> &collisions,
                             unsigned int moving_obj_type_id,
                             unsigned int still_obj_type_id,
						     scene_object obj1, scene_object obj2){        
    //TODO: need to handle different object differently? 
    //PLANE: check against the sides of the plane, just need plane normal and a plane vertex 
    //OTHERS: check against every vertices and compute distance -- no inside/outside check here? not needed?
    Eigen::VectorXd q = std::get<8>(obj1);
    Eigen::VectorXd q2 = std::get<8>(obj2);
    Eigen::MatrixXd sV = std::get<1>(obj2);
    Eigen::MatrixXi sF = std::get<2>(obj2);
	
    if (still_obj_type_id == 0) {
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
                collisions.push_back(std::make_pair(curr_v - dir * dist, vi));
            }
        }
    } else if (still_obj_type_id == -1)
    {
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

            double dist = (curr_v(1) - pos(1));
            if (dist <= 0.2) {
                collisions.push_back(std::make_pair(curr_v - dir * abs(dist), vi));
            }
        }
    } else
    {
        for (int vi = 0; vi < q.rows() / 3; vi++) {
            bool found = false;
            for (int fi = 0; fi < sF.rows(); fi++)
            {
                Eigen::Vector3d intersection_pt;
                found = point_triangle_intersection(q.segment<3>(vi * 3), q2, sF, intersection_pt);
            	if (found)
            	{
                    //std::cout << "\n" << intersection_pt << "\n";
                    collisions.push_back(std::make_pair(intersection_pt, vi));
            		break;
            	}
            }
        }
    }
}


void collision_detection(std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>> &collisions,
                         unsigned int moving_obj_geometry_idx,
                         unsigned int still_obj_geometry_idx,
                         scene_object obj1, scene_object obj2){       

    Eigen::VectorXd q = std::get<8>(obj1);
    Eigen::VectorXd qdot = std::get<9>(obj1);
    Eigen::VectorXd q2 = std::get<8>(obj2);
    Eigen::MatrixXd sV = std::get<1>(obj2);
    //std::cout<<std::get<0>(obj2)<<std::endl;
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
            double dist = (curr_v - pos).dot(dir);
            //double dist = abs((curr_v - pos).dot(dir));
            if (dist < 0.2) {
                //gather collision info to apply collision force
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

            Eigen::Vector3d dir;
            if (std::get<16>(obj2)(0) > 5)
            {
                dir = Eigen::Vector3d(-1, 0, 0);
            }
            else if (std::get<16>(obj2)(0) < -5)
            {
                dir = Eigen::Vector3d(1, 0, 0);
            }
            else if (std::get<16>(obj2)(2) < -5)
            {
                dir = Eigen::Vector3d(0, 0, 1);
            }
            else if (std::get<16>(obj2)(2) > 5)
            {
                dir = Eigen::Vector3d(0, 0, -1);
            }
            else
            {
                dir = Eigen::Vector3d(0, 1, 0);
            }
            dir = dir.normalized();
            double dist = (curr_v - pos).dot(dir);
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

void collision_detection_with_optimization(std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>>& collisions,
    unsigned int moving_obj_type_id,
    unsigned int still_obj_type_id,
    scene_object obj1, scene_object obj2) {
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
            double dist = (curr_v - pos).dot(dir);
            if (dist < 0.2) {
                std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int> collision_info;
                std::get<0>(collision_info) = curr_v - dir * dist;
                std::get<1>(collision_info) = -dir; // negative sign is added to make it the vertex normal
                std::get<2>(collision_info) = vi;
                std::get<3>(collision_info) = sF.row(0)(1);
                std::get<4>(collision_info) = still_obj_type_id;
                collisions.push_back(collision_info);
            }
        }
    }
    else if (std::get<0>(obj2) == -1)
    {
        for (int vi = 0; vi < q.rows() / 3; ++vi) {
            Eigen::Vector3d curr_v = q.segment<3>(3 * vi);
            //if it is plane 
            Eigen::Vector3d pos = sV.row(0);

            Eigen::Vector3d dir;
            if (std::get<16>(obj2)(0) > 5)
            {
                dir = Eigen::Vector3d(-1, 0, 0);
            }
            else if (std::get<16>(obj2)(0) < -5)
            {
                dir = Eigen::Vector3d(1, 0, 0);
            }
            else if (std::get<16>(obj2)(2) < -5)
            {
                dir = Eigen::Vector3d(0, 0, 1);
            }
            else if (std::get<16>(obj2)(2) > 5)
            {
                dir = Eigen::Vector3d(0, 0, -1);
            }
            else
            {
                dir = Eigen::Vector3d(0, 1, 0);
            }
            dir = dir.normalized();
            double dist = (curr_v - pos).dot(dir);
            if (dist <= 0.2) {

                std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int> collision_info;
                std::get<0>(collision_info) = curr_v - dir * dist;
                std::get<1>(collision_info) = -dir; // negative sign is added to make it the vertex normal
                std::get<2>(collision_info) = vi;
                std::get<3>(collision_info) = sF.row(0)(1);
                std::get<4>(collision_info) = still_obj_type_id;
                collisions.push_back(collision_info);
            }
        }
    }
    else
    {
        std::vector<std::vector<int>> VtF = std::get<17>(obj1);
        Eigen::MatrixXi F1 = std::get<2>(obj1);

    	std::vector<int> bb_list = std::get<19>(obj2);
        
        for (int vi = 0; vi < q.rows() / 3; vi++) {
            bool found_p2, found;
            found = false;
            int p2_index;
            found_p2 = false;
            Eigen::Vector3d p2, n1, p3, n2;
            double faces_count = (double)VtF.at(vi).size();
        	// calculate vertex normal
            for (int fi_vi = 0; fi_vi < faces_count; fi_vi++)
            {
                Eigen::Vector3d v0, v1, v2, n_of_vi;
                v0 = q.segment<3>(F1.row(VtF.at(vi).at(fi_vi))(0) * 3);
                v1 = q.segment<3>(F1.row(VtF.at(vi).at(fi_vi))(1) * 3);
                v2 = q.segment<3>(F1.row(VtF.at(vi).at(fi_vi))(2) * 3);
                n_of_vi = (v1 - v0).cross(v2 - v0);
                n1 = n1 + 1.0 / faces_count * n_of_vi;
            }
        	// getting the subset of faces
            std::vector<int> face_list_I_care_about;
            std::vector<int> tempA, tempB;
            std::unordered_map <int, Bounding_box> ht = std::get<18>(obj2);

            n1 = n1.normalized();
            for (int bi = 0; bi < bb_list.size(); bi++)
            {
                Bounding_box bb = ht.at(bb_list.at(bi));
            	if (bb.ray_box_intersection(q.segment<3>(3 * vi), n1) || bb.ray_box_intersection(q.segment<3>(3 * vi), -n1))
            	{
                    tempA = face_list_I_care_about;
                    bb.get_face_list(tempB);
                    face_list_I_care_about.reserve(tempA.size() + tempB.size());
                    face_list_I_care_about.insert(face_list_I_care_about.end(), tempA.begin(), tempA.end());
                    face_list_I_care_about.insert(face_list_I_care_about.end(), tempB.begin(), tempB.end());
            	}
            }

        	
            // see if there are any intersections between the vertex p1 on obj1 and faces on obj2
            for (int fici = 0; fici < face_list_I_care_about.size(); fici++)
            {
                int fi = face_list_I_care_about.at(fici);
                found = false;
                Eigen::Vector3d intersection_pt;
                double t;
                // obtain the vecters normal (using an average of the face normals around the vertex)
                if (found_p2) {
                    found = ray_triangle_intersect(p2, n1, t, q2, sF.row(fi));
                }
                else
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
                }
                else if (found && found_p2)
                {
                    p3 = q.segment<3>(vi * 3) + t * (n1);
                    if ((p2 - q.segment<3>(vi * 3)).norm() <= (p3 - q.segment<3>(vi * 3)).norm())
                    {
                        std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int> collision_info;
                        std::get<0>(collision_info) = p2;
                        std::get<1>(collision_info) = -n2;
                        std::get<2>(collision_info) = vi;
                        std::get<3>(collision_info) = p2_index;
                        std::get<4>(collision_info) = still_obj_type_id;
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
    } else if (std::get<0>(obj2) == 0)
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
        double dist = (com1 - pos).dot(dir);
    	return dist - rad1 <= 0.1;
    } else
    {

        Eigen::Vector3d dir;
        //std::cout << std::get<16>(obj2)(0) << "||" << std::get<16>(obj2)(1) << "||" << std::get<16>(obj2)(2) << "\n\n";
    	if (std::get<16>(obj2)(0) > 5)
    	{
            dir = Eigen::Vector3d(-1, 0, 0);
    	} else if (std::get<16>(obj2)(0) < -5)
    	{
            dir = Eigen::Vector3d(1, 0, 0);
    	} else if (std::get<16>(obj2)(2) < -5)
    	{
            dir = Eigen::Vector3d(0, 0, 1);
    	}
        else if (std::get<16>(obj2)(2) > 5)
        {
            dir = Eigen::Vector3d(0, 0, -1);
        } else
        {
            dir = Eigen::Vector3d(0, 1, 0);
        }
        //std::cout << dir(0) << "||" << dir(1) << "||" << dir(2) << "\n\n";
        com1 = std::get<16>(obj1);
        rad1 = std::get<15>(obj1);
        Eigen::MatrixXd sV = std::get<1>(obj2);
        Eigen::MatrixXi sF = std::get<2>(obj2);

        Eigen::Vector3d pos = sV.row(0).transpose();
        double dist = (com1 - pos).dot(dir);
        return dist - rad1 <= 0.1;
    }
}

void construct_spatial_hash_table(scene_object & obj)
{
    std::unordered_map <int, Bounding_box> ht;
    std::vector<int> occupied_hash_keys;
    Spatial_hash_fn hf;
    double cellsize = hf.cellsize;
    occupied_hash_keys.clear();
    ht.clear();
	
    Eigen::MatrixXi F = std::get<2>(obj);
    Eigen::VectorXd q = std::get<8>(obj);
    Eigen::Vector3d v0, v1, v2, max_cor, min_cor, max_cor_discreet, min_cor_discreet, temp_pt;
    for (int i = 0; i < F.rows(); i++)
    {
    	// obtain the maximum and minimum corner of the bound box
        v0 = q.segment<3>(F.row(i)(0) * 3);
        v1 = q.segment<3>(F.row(i)(1) * 3);
        v2 = q.segment<3>(F.row(i)(2) * 3);
        max_cor(0) = std::max(std::max(v0(0), v1(0)), v2(0));
        max_cor(1) = std::max(std::max(v0(1), v1(1)), v2(1));
        max_cor(2) = std::max(std::max(v0(2), v1(2)), v2(2));
        min_cor(0) = std::min(std::min(v0(0), v1(0)), v2(0));
        min_cor(1) = std::min(std::min(v0(1), v1(1)), v2(1));
        min_cor(2) = std::min(std::min(v0(2), v1(2)), v2(2));
        if (hf(min_cor) == hf(max_cor))
        {
        	// if the entire face lies in one boundbox
            if (ht.count(hf(min_cor)) == 0)
            {
            	// use the get_mapping function to get the discreet bbox corner
                hf.get_mapping(min_cor, temp_pt);
            	// create a bounding box and add it to the list of occupied boxes
                Bounding_box box = Bounding_box(temp_pt, temp_pt + Eigen::Vector3d::Ones() * cellsize);
                box.add_face(i);
                ht[hf(min_cor)] = box;
                occupied_hash_keys.push_back(hf(min_cor));
            } else
            {
                ht[hf(min_cor)].add_face(i);
            }
        } else
        {
            // if the entire face lies in multiple boundboxes
            hf.get_mapping(max_cor, max_cor_discreet);
            hf.get_mapping(min_cor, min_cor_discreet);
            int dx, dy, dz;
            double delta = cellsize / 10;
            dx = (int)round((max_cor_discreet(0) - min_cor_discreet(0)) / cellsize);
            dy = (int)round((max_cor_discreet(1) - min_cor_discreet(1)) / cellsize);
            dz = (int)round((max_cor_discreet(2) - min_cor_discreet(2)) / cellsize);
            for (int ix = 0; ix <= dx; ix++)
            {
                for (int iy = 0; iy <= dy; iy++)
                {
                    for (int iz = 0; iz <= dz; iz++)
                    {
                        temp_pt = min_cor_discreet;
                        temp_pt(0) += ix * cellsize;
                        temp_pt(1) += iy * cellsize;
                        temp_pt(2) += iz * cellsize;
                        if (ht.count(hf(temp_pt)) == 0)
                        {
                            Bounding_box box = Bounding_box(temp_pt, temp_pt + Eigen::Vector3d::Ones() * cellsize);
                        	box.add_face(i);
                            ht[hf(temp_pt)] = box;
                            occupied_hash_keys.push_back(hf(temp_pt));
                        }
                        else
                        {
                            ht[hf(temp_pt)].add_face(i);
                        }
                    	
                    }
                }
            }

        	
        }
    	
    }
    std::get<18>(obj) = ht;
    std::get<19>(obj) = occupied_hash_keys;
    return;
}
