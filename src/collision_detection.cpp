#include <collision_detection.h>
#include <iostream>
bool point_triangle_intersection(Eigen::Vector3d point, Eigen::VectorXd vertex_list, Eigen::RowVector3i face, Eigen::Vector3d & collision_pt)
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
//bool ray_triangle_intersect(
//    Eigen::Vector3d ray_origin, Eigen::Vector3d ray_direction, double& t, Eigen::VectorXd vertex_list, Eigen::RowVector3i face)
//{
//    double min_t = 0.0001;
//	
//    double local_t = 0;
//
//    
//	
//    Eigen::Vector3d s = vertex_list.row(face(0)).transpose() - ray_origin;
//    Eigen::Vector3d normal;
//    Eigen::Matrix3d C;
//    C.col(0) = vertex_list.row(face(0)).transpose() - vertex_list.row(face(1)).transpose();
//    C.col(1) = vertex_list.row(face(0)).transpose() - vertex_list.row(face(2)).transpose();
//    C.col(2) = ray_direction;
//    // pre calculating the variable before the computation. The equation from the text book is used.
//    double ei_hf, gf_di, dh_eg, ak_jb, jc_al, bl_kc;
//    double a = C(0, 0);
//    double b = C(1, 0);
//    double c = C(2, 0);
//    double d = C(0, 1);
//    double e = C(1, 1);
//    double f = C(2, 1);
//    double g = C(0, 2);
//    double h = C(1, 2);
//    double i = C(2, 2);
//    double j = s(0);
//    double k = s(1);
//    double l = s(2);
//    double beta, gamma;
//    ei_hf = e * i - h * f;
//    gf_di = g * f - d * i;
//    dh_eg = d * h - e * g;
//    double m = ei_hf * a + gf_di * b + dh_eg * c;
//    if (m == 0) {
//        return false;
//    }
//    ak_jb = a * k - j * b;
//    jc_al = j * c - a * l;
//    bl_kc = b * l - k * c;
//    // solve for t, gamma and beta sequentially so the algorithm can terminated earlier if needed
//    local_t = -(f * ak_jb + e * jc_al + d * bl_kc) / m;
//    if (local_t < min_t) {
//        return false;
//    }
//    gamma = (i * ak_jb + h * jc_al + g * bl_kc) / m;
//    if (gamma < 0 || gamma > 1) {
//        return false;
//    }
//    beta = (j * ei_hf + k * gf_di + l * dh_eg) / m;
//    if (beta < 0 || beta > 1 - gamma) {
//        return false;
//    }
//    // overwrite t with local value only if the ray actually intersect with the shape
//    t = local_t;
//    s = vertex_list.row(face(1)).transpose() - vertex_list.row(face(0)).transpose();
//    // computer normal and pick the normal that faces the eye (i.e. oposite to the ray from the eye)
//    //normal = s.cross3(vertex_list.row(face(2)).transpose() - vertex_list.row(face(0)).transpose());
//
//	return true;
//    ////////////////////////////////////////////////////////////////////////////
//}

void collision_detection(std::vector<std::pair<Eigen::Vector3d, unsigned int>> &collisions,
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
                    collisions.push_back(std::make_pair(intersection_pt, vi));
            		break;
            	}
            }
        }
    }
}
void compute_vertex_face_list(Eigen::MatrixXd V, Eigen::MatrixXd F, std::vector<std::vector<int>> &V2F)
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
