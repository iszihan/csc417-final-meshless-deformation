//#ifndef ASSIGNMENT_SETUP_H
//#define ASSIGNMENT_SETUP_H
//#include <tuple>
//#include <igl/readMESH.h>
//#include <igl/readOBJ.h>
//#include <igl/writeOBJ.h>
//#include <igl/readOFF.h>
//#include <igl/upsample.h>
//#include <read_tetgen.h>
//#include <igl/boundary_facets.h>
//#include <igl/volume.h>
////assignment files for implementing simulation and interaction
//#include <visualization.h>
//#include <init_state.h>
//#include <find_min_vertices.h>
//#include <fixed_point_constraints.h>
//#include <mass_matrix_particles.h>
//#include <meshless_implicit_euler.h>
//#include <dV_cloth_gravity_dq.h>
//#include <dV_spring_particle_particle_dq.h>
//#include <mutex>
//#include <collision_detection.h>
//
////variables for geometry
//typedef std::tuple<int, //moving or still
//Eigen::MatrixXd, //V 
//Eigen::MatrixXi, //F
//Eigen::MatrixXd, //V_skin
//Eigen::MatrixXi, //F_skin
//Eigen::SparseMatrixd, //N skinning matrix
//Eigen::SparseMatrixd, //M
//Eigen::Vector3d, //center of mass
//Eigen::VectorXd, //q
//Eigen::VectorXd, //qdot
//Eigen::SparseMatrixd, //P
//Eigen::VectorXd, //x0
//Eigen::VectorXd, //gravity
//std::vector<std::vector<int>>, //clusters
//Eigen::MatrixXd> scene_object;  //Q = V-center_of_mass
//
//std::string data_paths[3] = {"../data/cube.obj", 
//                             "../data/coarse_bunny2.obj",
//                             "../data/cube.obj"};
//std::vector<Eigen::VectorXd> force_list;
//
////material parameters
//double mass = 1.0;
//
////scratch memory for assembly
//Eigen::VectorXd qdot_tmp;
//Eigen::VectorXd tmp_force;
//Eigen::VectorXd q_tmp;
//Eigen::MatrixXd V_tmp;
//Eigen::MatrixXi F_tmp;
//Eigen::Vector3d com_tmp;
//int item_placement = -1;
//std::vector<std::vector<std::pair<Eigen::Vector3d, unsigned int>>> spring_points_list;
//std::vector<std::vector<std::pair<Eigen::Vector3d, unsigned int>>> collision_points_list;
////integration
//int method = 0;
//
////collision detection stuff
//bool collision_detection_on = false;
//bool simulation_pause = true;
//bool drop_key_pressed = false;
//int item_type = 0;
//
////selection spring
//double lspring = 0.1;
//double k_selected = 1e5;
//double k_collision = 1e5;
//
//inline void add_object_VF(std::vector<scene_object> &geometry, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool fixed, Eigen::Vector3i clusters)
//{
//    Eigen::VectorXd q;
//    Eigen::VectorXd qdot;
//    Eigen::VectorXd gravity;
//    Eigen::SparseMatrixd M;
//    Eigen::Vector3d center_of_mass;
//    Eigen::MatrixXd V_skin; //vertices of simulation mesh
//    Eigen::MatrixXi F_skin; //faces of simulation mesh
//    Eigen::SparseMatrixd N;
//    Eigen::SparseMatrixd P;
//    Eigen::VectorXd x0;
//    Eigen::MatrixXd Q;
//
//    init_state(q, qdot, V);
//    center_of_mass = V.colwise().mean();
//    Q = V.transpose().rowwise() - center_of_mass.transpose();
//
//    //skinning
//    V_skin = V;
//    F_skin = F;
//    N.resize(V.rows(), V.rows());
//    N.setIdentity();
//
//    //mass matrix
//    mass_matrix_particles(M, q, mass);
//    if (M.rows() == 0)
//    {
//        std::cout << "Mass matrix not implemented, exiting.\n";
//        exit(1);
//    }
//
//    //gravity vector
//    gravity.resize(q.rows(), 1);
//    dV_cloth_gravity_dq(gravity, M, Eigen::Vector3d(0, -9.8, 0));
//
//    if(fixed){
//        //fix to the floor
//        std::vector<unsigned int> fixed_point_indices;
//        find_min_vertices(fixed_point_indices, V, 0.001);
//        fixed_point_indices.push_back(0);
//        P.resize(q.rows(), q.rows());
//        P.setIdentity();
//        fixed_point_constraints(P, q.rows(), fixed_point_indices);
//        x0 = q - P.transpose() * P * q; //vector x0 contains position of all fixed nodes, zero for everything else
//        //correct M, q and qdot so they are the right size
//        // q = P * q;
//        // qdot = P * qdot;
//        // M = P * M * P.transpose();
//    }else{
//        //not fixed to the floor 
//        P.resize(q.rows(), q.rows());
//        P.setIdentity();
//        x0.resize(q.size());
//        x0.setZero();
//    }
//
//    // form clusters [x_clusters, y_cluster, z_clusters]
//    Eigen::Vector3d V_max = V.colwise().maxCoeff();
//    Eigen::Vector3d V_min = V.colwise().minCoeff();
//    double dx = (V_max(0)-V_min(0)) / clusters(0);
//    double dy = (V_max(1)-V_min(1)) / clusters(1);
//    double dz = (V_max(2)-V_min(2)) / clusters(2);
//    std::vector<std::vector<int>> vertex_clusters;
//    for(int ix=0; ix<clusters(0); ++ix){
//        double x_min = V_min(0) + dx * ix;
//        double x_max = V_min(0) + dx * (ix+1);
//        for(int iy=0; iy<clusters(1); ++iy){
//            double y_min = V_min(1) + dy * iy;
//            double y_max = V_min(1) + dy * (iy+1);
//            for(int iz=0; iz<clusters(2); ++iz){
//                double z_min = V_min(2) + dz * iz;
//                double z_max = V_min(2) + dz * (iz+1);
//                // std::cout<<"x range:"<<x_min<<"-"<<x_max<<std::endl;
//                // std::cout<<"y range:"<<y_min<<"-"<<y_max<<std::endl;
//                // std::cout<<"z range:"<<z_min<<"-"<<z_max<<std::endl;
//                std::vector<int> curr_cluster;
//                curr_cluster.clear();
//                double tol = 1e-3;
//                for(int iv=0; iv<V.rows();++iv){
//                    //put vertices to current cluster
//                    if(V.row(iv)(0) >= x_min - tol && V.row(iv)(0) <= x_max + tol
//                    && V.row(iv)(1) >= y_min - tol&& V.row(iv)(1) <= y_max + tol
//                    && V.row(iv)(2) >= z_min - tol&& V.row(iv)(2) <= z_max + tol){
//                        curr_cluster.push_back(iv);
//                    }
//                }
//                vertex_clusters.push_back(curr_cluster);
//            }
//        }
//    }
//    // for(int i=0; i<vertex_clusters.size();++i){
//    //     std::cout<<vertex_clusters.at(i).size()<<std::endl;
//    //     for(int j=0; j<vertex_clusters.at(i).size();++j){
//    //         std::cout<<V.row(vertex_clusters.at(i).at(j))<<std::endl;
//    //     }
//    // }
//
//    // append everything to the corresponding list
//    scene_object one_geometry;
//    std::get<0>(one_geometry) = 1;
//    std::get<1>(one_geometry) = V;
//    std::get<2>(one_geometry) = F;
//    std::get<3>(one_geometry) = V;
//    std::get<4>(one_geometry) = F;
//    std::get<5>(one_geometry) = N;
//    std::get<6>(one_geometry) = M;
//    std::get<7>(one_geometry) = center_of_mass;
//    std::get<8>(one_geometry) = q;
//    std::get<9>(one_geometry) = qdot;
//    std::get<10>(one_geometry) = P;
//    std::get<11>(one_geometry) = x0;
//    std::get<12>(one_geometry) = gravity;
//    std::get<13>(one_geometry) = vertex_clusters;
//    std::get<14>(one_geometry) = Q;
//	
//    geometry.push_back(one_geometry);
//    Visualize::add_object_to_scene(V, F, V, F, N, Eigen::RowVector3d(244, 165, 130) / 255.);
//}
//
//inline void add_object(std::vector<scene_object> &geometry, std::string file_path, Eigen::Vector3d position, bool fixed)
//{
//    Eigen::VectorXd q;
//    Eigen::VectorXd qdot;
//    Eigen::VectorXd gravity;
//    Eigen::SparseMatrixd M;
//    Eigen::Vector3d center_of_mass;
//    Eigen::MatrixXd V, V_skin; //vertices of simulation mesh
//    Eigen::MatrixXi F, F_skin; //faces of simulation mesh
//    Eigen::SparseMatrixd N;
//    Eigen::SparseMatrixd P;
//    Eigen::VectorXd x0;
//    Eigen::MatrixXd Q;
//
//    igl::readOBJ(file_path, V, F);
//
//    init_state(q, qdot, V);
//	for (int vi = 0; vi < q.size()/3; vi++)
//	{
//        q.segment<3>(vi * 3) += position;
//	}
//	
//    center_of_mass = V.colwise().mean();
//	
//    //add geometry to scene
//    V_skin = V;
//    F_skin = F;
//    N.resize(V.rows(), V.rows());
//    N.setIdentity();
//
//    //mass matrix
//    mass_matrix_particles(M, q, mass);
//    if (M.rows() == 0)
//    {
//        std::cout << "Mass matrix not implemented, exiting.\n";
//        exit(1);
//    }
//
//    //gravity vector
//    gravity.resize(q.rows(), 1);
//    dV_cloth_gravity_dq(gravity, M, Eigen::Vector3d(0, -9.8, 0));
//
//    if(fixed){
//        //fix to the floor
//        std::vector<unsigned int> fixed_point_indices;
//        find_min_vertices(fixed_point_indices, V, 0.001);
//        //fixed_point_indices.push_back(0);
//        P.resize(q.rows(), q.rows());
//        P.setIdentity();
//        fixed_point_constraints(P, q.rows(), fixed_point_indices);
//        x0 = q - P.transpose() * P * q; //vector x0 contains position of all fixed nodes, zero for everything else
//        //correct M, q and qdot so they are the right size
//        // q = P * q;
//        // qdot = P * qdot;
//        // M = P * M * P.transpose();
//    }else{
//        //not fixed to the floor 
//        P.resize(q.rows(), q.rows());
//        P.setIdentity();
//        x0.resize(q.size());
//        x0.setZero();
//    }
//
//    // append everything to the corresponding list
//
//    scene_object one_geometry;
//    std::get<0>(one_geometry) = 1;
//    std::get<1>(one_geometry) = V;
//    std::get<2>(one_geometry) = F;
//    std::get<3>(one_geometry) = V;
//    std::get<4>(one_geometry) = F;
//    std::get<5>(one_geometry) = N;
//    std::get<6>(one_geometry) = M;
//    std::get<7>(one_geometry) = center_of_mass;
//    std::get<8>(one_geometry) = q;
//    std::get<9>(one_geometry) = qdot;
//    std::get<10>(one_geometry) = P;
//    std::get<11>(one_geometry) = gravity;
//
//    geometry.push_back(one_geometry);
//}
//
//inline void simulate(std::vector<scene_object> &geometry, double dt, double t, std::mutex & mtx)
//{
//    //std::cout<<"inside simulate"<<std::endl;
//    force_list.clear();
//    //Interaction spring
//    if (!simulation_pause) {
//    	
//        spring_points_list.clear();
//        collision_points_list.clear();
//    	for (int i = 0; i < geometry.size(); i++)
//    	{
//            std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
//            std::vector<std::pair<Eigen::Vector3d, unsigned int>> collision_points_tmp;
//            spring_points_list.push_back(spring_points_tmp);
//            collision_points_list.push_back(collision_points_tmp);
//    	}
//
//        //add collision spring points with other geometry
//        Eigen::Vector6d dV_collide;
//        for (int mi = 0; mi < geometry.size(); ++mi) {
//            std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
//            std::vector<std::pair<Eigen::Vector3d, unsigned int>> collision_points_tmp;
//        	
//            scene_object moving_object = geometry.at(mi);
//        	// only check if the object is not static
//        	if (std::get<0>(moving_object) >= 1)
//        	{
//        		// check it against each other objects
//                for (unsigned int si = 0; si < geometry.size(); ++si) {
//                    // do not check collision against itself
//                    if (si != mi) {
//                    	// only consider if it is ball-parked to be touching
//                        scene_object collision_target = geometry.at(si);
//                    	//if (precomputation(moving_object, collision_target)){
//                        if (true) {
//                    		// only have collision detection with planes for now others can be implemented later
//                            if (std::get<0>(collision_target) == 0) {
//                                collision_detection(collision_points_list.at(mi), mi, si,
//                                    std::get<8>(moving_object), std::get<1>(collision_target), std::get<2>(collision_target));
//                            }
//                        }
//                    }
//                }
//        	}
//        }
//
//        Eigen::Vector3d mouse;
//        Eigen::Vector6d dV_mouse;
//        double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
//        for (int object_id = 0; object_id < geometry.size(); object_id++) {
//        	// only consider the movable objects to save compute
//            if (std::get<0>(geometry.at(object_id)) > 0) {
//                for (unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++)
//                {
//                    Eigen::VectorXd q = std::get<8>(geometry.at(object_id));
//                    Eigen::Vector3d p1 = (q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6);
//                    Eigen::Vector3d p2 = (q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]);
//                    spring_points_list.at(object_id).push_back(std::make_pair((q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi]));
//                    //TODO: add a dragging handle visualization
//                    //Visualize::scale_x(1, 2.0);
//                }
//            }
//        }
//        
//        for (int i = 0; i < geometry.size(); i++) {
//        //add collision spring points between moving geometries
//            scene_object current_object = geometry.at(i);
//            if (std::get<0>(current_object) > 0) {
//                auto force = [&](Eigen::VectorXd& f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2)
//                {
//                    Eigen::SparseMatrixd P = std::get<10>(current_object);
//                    std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points = spring_points_list.at(i);
//                    f = -std::get<12>(current_object);
//                    //dragging force
//                     //for (unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++)
//                     //{
//                     //    //dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose() * q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
//                     //    dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
//                     //    //std::cout<<"q1:"<<spring_points[pickedi].first<<std::endl;
//                     //    //std::cout<<"q2:"<<(q2 + x0).segment<3>(spring_points[pickedi].second)<<std::endl;
//                     //    f.segment<3>(3 * Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
//                     //    //std::cout<<"force:"<<std::endl;
//                     //    //std::cout<<dV_mouse.segment<3>(3)<<std::endl;
//                     //}
//
//                    //collision force
//                    for (unsigned int ci = 0; ci < collision_points_list.at(i).size(); ci++)
//                    {
//                        //dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose() * q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
//                        dV_spring_particle_particle_dq(dV_collide, collision_points_list.at(i)[ci].first, (q2 + std::get<11>(current_object)).segment<3>(3 * collision_points_list.at(i)[ci].second), 0.0, k_collision);
//                        //std::cout << "q1:" << collision_points[ci].first << std::endl;
//                        //std::cout << "q2:" << (q2 + x0).segment<3>(3 * collision_points[ci].second) << std::endl;
//                        f.segment<3>(3 * collision_points_list.at(i)[ci].second) += dV_collide.segment<3>(3);
//                        //std::cout << "added force:" << std::endl;
//                        //std::cout << dV_collide.segment<3>(3) << std::endl;
//                        //std::cout << "current force:" << std::endl;
//                        //std::cout << f.segment<3>(3 * collision_points[ci].second) << std::endl;
//                    }
//
//                    f = P * f;
//                };
//
//                q_tmp = std::get<8>(current_object);
//                qdot_tmp = std::get<9>(current_object);
//                V_tmp = std::get<1>(current_object);
//                F_tmp = std::get<2>(current_object);
//                com_tmp = std::get<7>(current_object);
//
//                
//                meshless_implicit_euler(q_tmp, qdot_tmp, dt, mass, V_tmp, com_tmp, force, tmp_force);
//                std::get<8>(geometry.at(i)) = q_tmp;
//                std::get<9>(geometry.at(i)) = qdot_tmp;
//                //std::get<13>(geometry.at(i)) = com_tmpt;
//            }
//		}
//
//        if (item_placement >= 0)
//        {
//            mtx.lock();
//            Eigen::Vector3d pos = Visualize::mouse_world();
//            //pos(1) = std::min(pos(1), 3.0);
//            add_object(geometry, data_paths[item_placement], pos, false);
//            item_placement = -1;
//            mtx.unlock();
//        }
//    }
//}
//
//inline void draw(std::vector<scene_object> geometry, double t)
//{
//    //update vertex positions using simulation
//    for (int i = 0; i < geometry.size(); i++) {
//        scene_object current_object = geometry.at(i);
//        if (std::get<0>(current_object) > 0) {
//		   Visualize::update_vertex_positions(i, std::get<10>(current_object).transpose() * std::get<8>(current_object) + std::get<11>(current_object));
//		}
//    }
//}
//
//bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers)
//{
//    if (key == 'N')
//    {
//        std::cout << "toggle integrators \n";
//        //fully_implicit = !fully_implicit;
//    }
//    if (key == 'L')
//    {
//        std::cout << "toggle deformation to linear \n";
//        method = 1;
//    }
//    if (key == 'Q')
//    {
//        std::cout << "toggle deformation to quadratic \n";
//        method = 2;
//    }if (key == 'R')
//    {
//        std::cout << "toggle deformation to rigid only \n";
//        method = 0;
//    }
//    if (key == 'C')
//    {
//        collision_detection_on = !collision_detection_on;
//        Visualize::set_visible(1, collision_detection_on);
//    }
//    //TODO: 
//    // - toggle key to drop objects from sky
//    if (key == 'Q'){
//        std::exit(1);
//    }
//    if (key == 'P'){
//        simulation_pause = !simulation_pause;
//    }
//
//	if ((key == '1' || key == '2' || key == '3') && item_placement == -1)
//	{
//		if (key == '1'){
//            item_placement = 0;
//        } else if (key == '2') {
//            item_placement = 1;
//        } else if (key == '3') {
//            item_placement = 2;
//        }
//        std::cout << "\n gonna be dropping item " << item_placement << "\n";
//	}
//	
//    return false;
//}
//inline void add_plane(Eigen::Vector3d floor_normal, Eigen::Vector3d floor_pos, std::vector<scene_object> & scene_object_list)
//{
//    Eigen::MatrixXd V_floor;
//    Eigen::MatrixXi F_floor;
//    Eigen::SparseMatrixd N;
//    igl::readOBJ("../data/plane.obj", V_floor, F_floor);
//    V_floor *= 10.0; //make it bigger
//    N.resize(V_floor.rows(), V_floor.rows());
//    N.setIdentity();
//    //rotate plane
//    Eigen::Vector3d n0;
//    n0 << 0, 1, 0;
//    floor_normal.normalize();
//    float angle = std::acos(n0.dot(floor_normal)) + 0.15;
//    Eigen::Vector3d axis = n0.cross(floor_normal);
//    Eigen::Matrix3d floor_R = Eigen::AngleAxisd(angle, axis).matrix();
//    //translate plane
//    for (unsigned int iv = 0; iv < V_floor.rows(); ++iv) {
//        Eigen::Vector3d rotated = floor_R * V_floor.row(iv).transpose();
//        V_floor.row(iv) = (rotated + floor_pos + Eigen::Vector3d(0., -0.01, 0.)).transpose();
//    }
//    
//    Visualize::add_object_to_scene(V_floor, F_floor, V_floor, F_floor, N, Eigen::RowVector3d(64, 165, 130) / 255.);
//    scene_object one_geometry;
//    Eigen::SparseMatrixd M;
//    std::get<0>(one_geometry) = 0;
//    std::get<1>(one_geometry) = V_floor;
//    std::get<2>(one_geometry) = F_floor;
//    std::get<3>(one_geometry) = V_floor;
//    std::get<4>(one_geometry) = F_floor;
//    std::get<5>(one_geometry) = N;
//    std::get<6>(one_geometry) = M;
//    std::get<7>(one_geometry) = Eigen::Vector3d::Zero();
//    std::get<8>(one_geometry) = Eigen::Vector3d::Zero();
//    std::get<9>(one_geometry) = Eigen::Vector3d::Zero();
//    scene_object_list.push_back(one_geometry);
//}
//inline void assignment_setup(int argc, char **argv, std::vector<Eigen::VectorXd> &q_list, std::vector<Eigen::VectorXd>&qdot_list, std::vector<scene_object> &geometry)
//{
//	
//    add_object()
//
//	
//    //add floor
//    Eigen::Vector3d floor_normal;
//    Eigen::Vector3d floor_pos;
//    //Eigen::SparseMatrixd N;
//    floor_normal << 0.0, 1.0, 0.0;
//    floor_pos << 0.0, -1.0, 0.0;    
//    add_plane(floor_normal, floor_pos, geometry);
//	
//    //Eigen::Vector3d floor_normal2;
//    //Eigen::Vector3d floor_pos2;
//    //floor_normal2 << 1, 0.7, 0.;
//    //floor_pos2 << 0., -3.0, 0.;
//    //add_plane(floor_normal2, floor_pos, geometry);
//    
//    Visualize::viewer().callback_key_down = key_down_callback;
//}
//#endif



#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H
#include <tuple>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
#include <igl/upsample.h>
#include <read_tetgen.h>
#include <igl/boundary_facets.h>
#include <igl/volume.h>
//assignment files for implementing simulation and interaction
#include <visualization.h>
#include <init_state.h>
#include <find_min_vertices.h>
#include <fixed_point_constraints.h>
#include <mass_matrix_particles.h>
#include <meshless_implicit_euler.h>
#include <dV_cloth_gravity_dq.h>
#include <dV_spring_particle_particle_dq.h>
#include <mutex>
#include <collision_detection.h>

//variables for geometry
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
				   Eigen::Vector3d                 // com that moves with time
                   > scene_object;
std::string data_paths[3] = {"../data/cube.obj",
                             "../data/coarse_bunny2.obj",
                             "../data/cube.obj"};

//material parameters
double mass = 1.0;

//scratch memory for assembly
Eigen::VectorXd qdot_tmp;
Eigen::VectorXd tmp_force;
Eigen::VectorXd q_tmp;
Eigen::MatrixXd V_tmp;
Eigen::MatrixXi F_tmp;
Eigen::Vector3d com_tmp;
Eigen::SparseMatrixd P_tmp;
Eigen::VectorXd x0_tmp;
Eigen::MatrixXd Q_tmp;
int item_placement = -1;
std::vector<std::vector<std::pair<Eigen::Vector3d, unsigned int>>> spring_points_list;
std::vector<std::vector<std::pair<Eigen::Vector3d, unsigned int>>> collision_points_list;
//integration
int method = 2;

//collision detection stuff
bool collision_detection_on = false;
bool simulation_pause = true;
bool drop_key_pressed = false;
int item_type = 0;

//selection spring
double lspring = 0.1;
double k_selected = 1e5;
double k_collision = 1e5;

inline void add_object_VF(std::vector<scene_object> &geometry, Eigen::MatrixXd &V, Eigen::MatrixXi &F, bool fixed, Eigen::Vector3i clusters)
{
    Eigen::VectorXd q;
    Eigen::VectorXd qdot;
    Eigen::VectorXd gravity;
    Eigen::SparseMatrixd M;
    Eigen::Vector3d center_of_mass;
    Eigen::MatrixXd V_skin; //vertices of simulation mesh
    Eigen::MatrixXi F_skin; //faces of simulation mesh
    Eigen::SparseMatrixd N;
    Eigen::SparseMatrixd P;
    Eigen::VectorXd x0;
    Eigen::MatrixXd Q;

    init_state(q, qdot, V);
    center_of_mass = V.colwise().mean();
    Q = V.rowwise() - center_of_mass.transpose();

    //skinning
    V_skin = V;
    F_skin = F;
    N.resize(V.rows(), V.rows());
    N.setIdentity();

    //mass matrix
    mass_matrix_particles(M, q, mass);
    if (M.rows() == 0)
    {
        std::cout << "Mass matrix not implemented, exiting.\n";
        exit(1);
    }

    //gravity vector
    gravity.resize(q.rows(), 1);
    dV_cloth_gravity_dq(gravity, M, Eigen::Vector3d(0, -9.8, 0));

    if (fixed)
    {
        //fix to the floor
        std::vector<unsigned int> fixed_point_indices;
        find_min_vertices(fixed_point_indices, V, 0.001);
        fixed_point_indices.push_back(0);
        P.resize(q.rows(), q.rows());
        P.setIdentity();
        fixed_point_constraints(P, q.rows(), fixed_point_indices);
        x0 = q - P.transpose() * P * q; //vector x0 contains position of all fixed nodes, zero for everything else
        //correct M, q and qdot so they are the right size
        // q = P * q;
        // qdot = P * qdot;
        // M = P * M * P.transpose();
    }
    else
    {
        //not fixed to the floor
        P.resize(q.rows(), q.rows());
        P.setIdentity();
        x0.resize(q.size());
        x0.setZero();
    }

    // form clusters [x_clusters, y_cluster, z_clusters]
    Eigen::Vector3d V_max = V.colwise().maxCoeff();
    Eigen::Vector3d V_min = V.colwise().minCoeff();
    double dx = (V_max(0) - V_min(0)) / clusters(0);
    double dy = (V_max(1) - V_min(1)) / clusters(1);
    double dz = (V_max(2) - V_min(2)) / clusters(2);
    std::vector<std::vector<int>> vertex_clusters;
    for (int ix = 0; ix < clusters(0); ++ix)
    {
        double x_min = V_min(0) + dx * ix;
        double x_max = V_min(0) + dx * (ix + 1);
        for (int iy = 0; iy < clusters(1); ++iy)
        {
            double y_min = V_min(1) + dy * iy;
            double y_max = V_min(1) + dy * (iy + 1);
            for (int iz = 0; iz < clusters(2); ++iz)
            {
                double z_min = V_min(2) + dz * iz;
                double z_max = V_min(2) + dz * (iz + 1);
                // std::cout<<"x range:"<<x_min<<"-"<<x_max<<std::endl;
                // std::cout<<"y range:"<<y_min<<"-"<<y_max<<std::endl;
                // std::cout<<"z range:"<<z_min<<"-"<<z_max<<std::endl;
                std::vector<int> curr_cluster;
                curr_cluster.clear();
                double tol = 1e-3;
                for (int iv = 0; iv < V.rows(); ++iv)
                {
                    //put vertices to current cluster
                    if (V.row(iv)(0) >= x_min - tol && V.row(iv)(0) <= x_max + tol && V.row(iv)(1) >= y_min - tol && V.row(iv)(1) <= y_max + tol && V.row(iv)(2) >= z_min - tol && V.row(iv)(2) <= z_max + tol)
                    {
                        curr_cluster.push_back(iv);
                    }
                }
                vertex_clusters.push_back(curr_cluster);
            }
        }
    }
    // for(int i=0; i<vertex_clusters.size();++i){
    //     std::cout<<vertex_clusters.at(i).size()<<std::endl;
    //     for(int j=0; j<vertex_clusters.at(i).size();++j){
    //         std::cout<<V.row(vertex_clusters.at(i).at(j))<<std::endl;
    //     }
    // }

    // append everything to the corresponding list
    scene_object one_geometry;
    std::get<0>(one_geometry) = 1;
    std::get<1>(one_geometry) = V;
    std::get<2>(one_geometry) = F;
    std::get<3>(one_geometry) = V;
    std::get<4>(one_geometry) = F;
    std::get<5>(one_geometry) = N;
    std::get<6>(one_geometry) = M;
    std::get<7>(one_geometry) = center_of_mass;
    std::get<8>(one_geometry) = q;
    std::get<9>(one_geometry) = qdot;
    std::get<10>(one_geometry) = P;
    std::get<11>(one_geometry) = x0;
    std::get<12>(one_geometry) = gravity;
    std::get<13>(one_geometry) = vertex_clusters;
    std::get<14>(one_geometry) = Q;
    geometry.push_back(one_geometry);
    Visualize::add_object_to_scene(V, F, V, F, N, Eigen::RowVector3d(244, 165, 130) / 255.);
}

inline void add_object(std::vector<scene_object>& geometry, std::string file_path, Eigen::Vector3d position, bool fixed)
{
    Eigen::VectorXd q;
    Eigen::VectorXd qdot;
    Eigen::VectorXd gravity;
    Eigen::SparseMatrixd M;
    Eigen::Vector3d center_of_mass;
    Eigen::MatrixXd V, V_skin; //vertices of simulation mesh
    Eigen::MatrixXi F, F_skin; //faces of simulation mesh
    Eigen::SparseMatrixd N;
    Eigen::SparseMatrixd P;
    Eigen::VectorXd x0;
    Eigen::MatrixXd Q;

    igl::readOBJ(file_path, V, F);

    init_state(q, qdot, V);
    for (int vi = 0; vi < q.size() / 3; vi++)
    {
        q.segment<3>(vi * 3) += position;
    }
    center_of_mass = V.colwise().mean();
    Q = V.rowwise() - center_of_mass.transpose();

    //skinning
    V_skin = V;
    F_skin = F;
    N.resize(V.rows(), V.rows());
    N.setIdentity();

    //mass matrix
    mass_matrix_particles(M, q, mass);
    if (M.rows() == 0)
    {
        std::cout << "Mass matrix not implemented, exiting.\n";
        exit(1);
    }

    //gravity vector
    gravity.resize(q.rows(), 1);
    dV_cloth_gravity_dq(gravity, M, Eigen::Vector3d(0, -9.8, 0));

    if (fixed)
    {
        //fix to the floor
        std::vector<unsigned int> fixed_point_indices;
        find_min_vertices(fixed_point_indices, V, 0.001);
        //fixed_point_indices.push_back(0);
        P.resize(q.rows(), q.rows());
        P.setIdentity();
        fixed_point_constraints(P, q.rows(), fixed_point_indices);
        x0 = q - P.transpose() * P * q; //vector x0 contains position of all fixed nodes, zero for everything else
        //correct M, q and qdot so they are the right size
        // q = P * q;
        // qdot = P * qdot;
        // M = P * M * P.transpose();
    }
    else
    {
        //not fixed to the floor
        P.resize(q.rows(), q.rows());
        P.setIdentity();
        x0.resize(q.size());
        x0.setZero();
    }

    //form a single cluster
    std::vector<std::vector<int>> vertex_clusters;
    std::vector<int> all_vertices(V.rows());
    std::iota(all_vertices.begin(), all_vertices.end(), V.rows());
    vertex_clusters.push_back(all_vertices);

    // append everything to the corresponding list
    scene_object one_geometry;
    std::get<0>(one_geometry) = 1;
    std::get<1>(one_geometry) = V;
    std::get<2>(one_geometry) = F;
    std::get<3>(one_geometry) = V;
    std::get<4>(one_geometry) = F;
    std::get<5>(one_geometry) = N;
    std::get<6>(one_geometry) = M;
    std::get<7>(one_geometry) = center_of_mass;
    std::get<8>(one_geometry) = q;
    std::get<9>(one_geometry) = qdot;
    std::get<10>(one_geometry) = P;
    std::get<11>(one_geometry) = x0;
    std::get<12>(one_geometry) = gravity;
    std::get<13>(one_geometry) = vertex_clusters;
    std::get<14>(one_geometry) = Q;
    double radius = ((V.rowwise() - center_of_mass.transpose()).rowwise().norm()).maxCoeff() * 1.5;
    std::get<15>(one_geometry) = radius;
    std::get<16>(one_geometry) = center_of_mass;

    geometry.push_back(one_geometry);
    Visualize::add_object_to_scene(V, F, V, F, N, Eigen::RowVector3d(244, 165, 130) / 255.);
}
inline void add_plane(Eigen::Vector3d floor_normal, Eigen::Vector3d floor_pos, std::vector<scene_object>& geometry)
{
    Eigen::MatrixXd V_floor;
    Eigen::MatrixXi F_floor;
    Eigen::SparseMatrixd N;
    igl::readOBJ("../data/plane.obj", V_floor, F_floor);
    //make it bigger
    V_floor *= 10.0;
    //skinning
    N.resize(V_floor.rows(), V_floor.rows());
    N.setIdentity();
    //rotate plane
    Eigen::Vector3d n0;
    n0 << 0, 1, 0;
    floor_normal.normalize();
    float angle = std::acos(n0.dot(floor_normal)) + 0.15;
    Eigen::Vector3d axis = n0.cross(floor_normal);
    Eigen::Matrix3d floor_R = Eigen::AngleAxisd(angle, axis).matrix();
    //translate plane
    for (unsigned int iv = 0; iv < V_floor.rows(); ++iv)
    {
        Eigen::Vector3d rotated = floor_R * V_floor.row(iv).transpose();
        V_floor.row(iv) = (rotated + floor_pos + Eigen::Vector3d(0., -0.01, 0.)).transpose();
    }

    Visualize::add_object_to_scene(V_floor, F_floor, V_floor, F_floor, N, Eigen::RowVector3d(64, 165, 130) / 255.);
    scene_object one_geometry;
    Eigen::SparseMatrixd M;
    std::get<0>(one_geometry) = 0;
    std::get<1>(one_geometry) = V_floor;
    std::get<2>(one_geometry) = F_floor;
    std::get<3>(one_geometry) = V_floor;
    std::get<4>(one_geometry) = F_floor;
    std::get<5>(one_geometry) = N;
    std::get<6>(one_geometry) = M;
    std::get<7>(one_geometry) = Eigen::Vector3d::Zero();
    std::get<8>(one_geometry) = Eigen::Vector3d::Zero();
    std::get<9>(one_geometry) = Eigen::Vector3d::Zero();
    geometry.push_back(one_geometry);
}
inline void add_floor(Eigen::Vector3d floor_normal, Eigen::Vector3d floor_pos, std::vector<scene_object>& geometry)
{
    Eigen::MatrixXd V_floor;
    Eigen::MatrixXi F_floor;
    Eigen::SparseMatrixd N;
    igl::readOBJ("../data/plane.obj", V_floor, F_floor);
    //make it bigger
    V_floor *= 10.0;
    //skinning
    N.resize(V_floor.rows(), V_floor.rows());
    N.setIdentity();
    //rotate plane
    Eigen::Vector3d n0;
    n0 << 0, 1, 0;
    floor_normal.normalize();
    float angle = std::acos(n0.dot(floor_normal)) + 0.15;
    Eigen::Vector3d axis = n0.cross(floor_normal);
    Eigen::Matrix3d floor_R = Eigen::AngleAxisd(angle, axis).matrix();
    //translate plane
    for (unsigned int iv = 0; iv < V_floor.rows(); ++iv)
    {
        Eigen::Vector3d rotated = floor_R * V_floor.row(iv).transpose();
        V_floor.row(iv) = (rotated + floor_pos + Eigen::Vector3d(0., -0.01, 0.)).transpose();
    }

    Visualize::add_object_to_scene(V_floor, F_floor, V_floor, F_floor, N, Eigen::RowVector3d(64, 165, 130) / 255.);
    scene_object one_geometry;
    Eigen::SparseMatrixd M;
    std::get<0>(one_geometry) = -1;
    std::get<1>(one_geometry) = V_floor;
    std::get<2>(one_geometry) = F_floor;
    std::get<3>(one_geometry) = V_floor;
    std::get<4>(one_geometry) = F_floor;
    std::get<5>(one_geometry) = N;
    std::get<6>(one_geometry) = M;
    std::get<7>(one_geometry) = Eigen::Vector3d::Zero();
    std::get<8>(one_geometry) = Eigen::Vector3d::Zero();
    std::get<9>(one_geometry) = Eigen::Vector3d::Zero();
    geometry.push_back(one_geometry);
}

inline void simulate(std::vector<scene_object> &geometry, double dt, double t, std::mutex &mtx)
{
    std::cout<<"inside simulate"<<std::endl;
    //Interaction spring
    if (!simulation_pause)
    {

        spring_points_list.clear();
        collision_points_list.clear();
        for (int i = 0; i < geometry.size(); i++)
        {
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> collision_points_tmp;
            spring_points_list.push_back(spring_points_tmp);
            collision_points_list.push_back(collision_points_tmp);
        }

        //add collision spring points with other geometry
        Eigen::Vector6d dV_collide;
        for (int mi = 0; mi < geometry.size(); ++mi)
        {
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> collision_points_tmp;
            scene_object moving_object = geometry.at(mi);
            // only check if the object is not static
            if (std::get<0>(moving_object) >= 1)
            {
                for (unsigned int si = 0; si < geometry.size(); ++si)
                {
                    // do not check collision against itself
                    if (si != mi)
                    {
                        scene_object collision_target = geometry.at(si);
                    	// only do collision compute if they overlapps
                        if (precomputation(moving_object, collision_target)) {
                            // only have collision detection with planes for now others can be implemented later
                            //std::cout << "checking collision" << std::endl;
                            collision_detection(collision_points_list.at(mi), std::get<0>(moving_object), std::get<0>(collision_target),
                                moving_object, collision_target);
                        }
                    }
                }
            }

        }
        Eigen::Vector3d mouse;
        Eigen::Vector6d dV_mouse;
        double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
        for (int object_id = 0; object_id < geometry.size(); object_id++)
        {
            // only consider the movable objects to save compute
            if (std::get<0>(geometry.at(object_id)) > 0)
            {
                for (unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++)
                {
                    Eigen::VectorXd q = std::get<8>(geometry.at(object_id));
                    Eigen::Vector3d p1 = (q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6);
                    Eigen::Vector3d p2 = (q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]);
                    spring_points_list.at(object_id).push_back(std::make_pair((q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi]));
                    //TODO: add a dragging handle visualization
                    //Visualize::scale_x(1, 2.0);
                }
            }
        }

        for (int i = 0; i < geometry.size(); i++)
        {
            //add collision spring points between moving geometries
            scene_object current_object = geometry.at(i);
            if (std::get<0>(current_object) > 0)
            {
                auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2)
                {
                    Eigen::SparseMatrixd P = std::get<10>(current_object);
                    std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points = spring_points_list.at(i);
                    f = -std::get<12>(current_object);
                    //dragging force
                    //for (unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++)
                    //{
                    //    //dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose() * q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                    //    dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                    //    //std::cout<<"q1:"<<spring_points[pickedi].first<<std::endl;
                    //    //std::cout<<"q2:"<<(q2 + x0).segment<3>(spring_points[pickedi].second)<<std::endl;
                    //    f.segment<3>(3 * Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
                    //    //std::cout<<"force:"<<std::endl;
                    //    //std::cout<<dV_mouse.segment<3>(3)<<std::endl;
                    //}

                    //collision force
                    for (unsigned int ci = 0; ci < collision_points_list.at(i).size(); ci++)
                    {
                        //dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose() * q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                        dV_spring_particle_particle_dq(dV_collide, collision_points_list.at(i)[ci].first, (q2 + std::get<11>(current_object)).segment<3>(3 * collision_points_list.at(i)[ci].second), 0.0, k_collision);
                        //std::cout << "q1:" << collision_points[ci].first << std::endl;
                        //std::cout << "q2:" << (q2 + x0).segment<3>(3 * collision_points[ci].second) << std::endl;
                        f.segment<3>(3 * collision_points_list.at(i)[ci].second) += dV_collide.segment<3>(3);
                        //std::cout << "added force:" << std::endl;
                        //std::cout << dV_collide.segment<3>(3) << std::endl;
                        //std::cout << "current force:" << std::endl;
                        //std::cout << f.segment<3>(3 * collision_points[ci].second) << std::endl;
                    }

                    f = P * f;
                };
                Eigen::Vector3d comt;
            	
                q_tmp = std::get<8>(current_object);
                qdot_tmp = std::get<9>(current_object);
                V_tmp = std::get<1>(current_object);
                F_tmp = std::get<2>(current_object);
                com_tmp = std::get<7>(current_object);
                P_tmp = std::get<10>(current_object);
                x0_tmp = std::get<11>(current_object);
                Q_tmp = std::get<14>(current_object);
                meshless_implicit_euler(q_tmp, qdot_tmp, dt, mass, P_tmp, x0_tmp, Q_tmp, std::get<13>(current_object), method, force, tmp_force, comt);
                std::get<8>(geometry.at(i)) = q_tmp;
                std::get<9>(geometry.at(i)) = qdot_tmp;
                std::get<16>(geometry.at(i)) = comt;
            }
        }

        if (item_placement >= 0)
        {
            mtx.lock();
            Eigen::Vector3d pos = Visualize::mouse_world();
            pos(1) = 3.0;
            pos(2) = 0.0;
            //pos(1) = std::min(pos(1), 3.0);
            add_object(geometry, data_paths[item_placement], pos, false);
            item_placement = -1;
            mtx.unlock();
        }
    }
}

inline void draw(std::vector<scene_object> geometry, double t)
{
    //update vertex positions using simulation
    for (int i = 0; i < geometry.size(); i++)
    {
        scene_object current_object = geometry.at(i);
        if (std::get<0>(current_object) > 0)
        {
            Visualize::update_vertex_positions(i, std::get<10>(current_object).transpose() * std::get<8>(current_object) + std::get<11>(current_object));
        }
    }
}

bool key_down_callback(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
    if (key == 'N')
    {
        std::cout << "toggle integrators \n";
        //fully_implicit = !fully_implicit;
    }
    if (key == 'L')
    {
        std::cout << "toggle deformation to linear \n";
        method = 1;
    }
    if (key == 'Q')
    {
        std::cout << "toggle deformation to quadratic \n";
        method = 2;
    }
    if (key == 'R')
    {
        std::cout << "toggle deformation to rigid only \n";
        method = 0;
    }
    if (key == 'C')
    {
        collision_detection_on = !collision_detection_on;
        Visualize::set_visible(1, collision_detection_on);
    }
    //TODO:
    // - toggle key to drop objects from sky
    if (key == 'E')
    {
        std::exit(1);
    }
    if (key == 'P')
    {
        simulation_pause = !simulation_pause;
    }

    if ((key == '1' || key == '2' || key == '3') && item_placement == -1)
    {
        if (key == '1')
        {
            item_placement = 0;
        }
        else if (key == '2')
        {
            item_placement = 1;
        }
        else if (key == '3')
        {
            item_placement = 2;
        }
        std::cout << "\n gonna be dropping item " << item_placement << "\n";
    }

    return false;
}

inline void simulate_clustering(std::vector<scene_object> &geometry, double dt, double t)
{
    //Interaction spring
    if (!simulation_pause)
    {
        //collect dragging points
        spring_points_list.clear();
        for (int i = 0; i < geometry.size(); i++)
        {
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
            spring_points_list.push_back(spring_points_tmp);
        }

        Eigen::Vector3d mouse;
        Eigen::Vector6d dV_mouse;
        double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
        for (int object_id = 0; object_id < geometry.size(); object_id++)
        {
            // only consider the movable objects to save compute
            if (std::get<0>(geometry.at(object_id)) > 0)
            {
                for (unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++)
                {
                    Eigen::VectorXd q = std::get<8>(geometry.at(object_id));
                    Eigen::Vector3d p1 = (q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6);
                    Eigen::Vector3d p2 = (q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]);
                    spring_points_list.at(object_id).push_back(std::make_pair((q + std::get<11>(geometry.at(object_id))).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi]));

                    //TODO: add a dragging handle visualization
                }
            }
        }

        //apply forces
        for (int i = 0; i < geometry.size(); i++)
        {
            scene_object current_object = geometry.at(i);
            if (std::get<0>(current_object) > 0)
            {
                auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2)
                {
                    //gravity
                    //f = -std::get<12>(current_object);
                    Eigen::SparseMatrixd P = std::get<10>(current_object);
                    std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points = spring_points_list.at(i);
                    for (unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++)
                    {
                        dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (q2 + std::get<11>(current_object)).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                        //std::cout<<"q1:"<<spring_points[pickedi].first<<std::endl;
                        //std::cout<<"q2:"<<(q2 + x0).segment<3>(spring_points[pickedi].second)<<std::endl;
                        f.segment<3>(3 * Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
                        //std::cout<<"force:"<<std::endl;
                        //std::cout<<dV_mouse.segment<3>(3)<<std::endl;
                    }
                    f = P * f;
                };

                q_tmp = std::get<8>(current_object);
                qdot_tmp = std::get<9>(current_object);
                V_tmp = std::get<1>(current_object);
                F_tmp = std::get<2>(current_object);
                com_tmp = std::get<7>(current_object);
                P_tmp = std::get<10>(current_object);
                x0_tmp = std::get<11>(current_object);
                Q_tmp = std::get<14>(current_object);
                Eigen::Vector3d comt;
                meshless_implicit_euler(q_tmp, qdot_tmp, dt, mass, P_tmp, x0_tmp, Q_tmp, std::get<13>(current_object), method, force, tmp_force, comt);
                std::get<8>(geometry.at(i)) = q_tmp;
                std::get<9>(geometry.at(i)) = qdot_tmp;
            }
        }
    }
}

inline void assignment_setup(int argc, char **argv, std::vector<scene_object> &geometry)
{
<<<<<<< HEAD

    Eigen::Vector3d origin;
    origin << 0.0, 5, 0.0;
    add_object(geometry, "../data/coarse_bunny2.obj", origin, false);

    origin << 0.0, 0, 0.0;
    add_object(geometry, "../data/coarse_bunny2.obj", origin, false);
=======
    // // load setup scene
    // Eigen::Vector3d origin;
    // origin << 0.0, 0.0, 0.0;
    // add_object(geometry, "../data/cube.obj", origin, false);

    //load in cube
    Eigen::MatrixXd V, SV;
    Eigen::MatrixXi F, SF;
    int subdiv = 2;
    igl::readOBJ("../data/cube.obj", V, F);
    //subdivide by 2
    for (int i = 0; i < subdiv; ++i)
    {
        igl::upsample(V, F, SV, SF);
        V = SV;
        F = SF;
    }
    Eigen::Vector3i cluster_size;
    cluster_size << 1, 1, 1;
    add_object_VF(geometry, SV, SF, false, cluster_size);

>>>>>>> ca1c9d8c1ad74c70a4cca7a2297c4639fe0315f6
    Eigen::Vector3d floor_normal;
    Eigen::Vector3d floor_pos;
    Eigen::SparseMatrixd N;
    floor_normal << 0.0, 1.0, 0.0;
    floor_pos << 0.0, -4.0, 0.0;
    //add_plane(floor_normal, floor_pos, geometry);
    add_floor(floor_normal, floor_pos, geometry);
    Visualize::viewer().callback_key_down = key_down_callback;
    //simulation_pause = false;
    std::cout << "finished set up" << std::endl;
    std::cout << "there is a total of " << geometry.size() << " pieces of geometry. \n";
}


inline void clustering_setup(int argc, char **argv, std::vector<scene_object> &geometry)
{
    //load in cube
    Eigen::MatrixXd V, SV;
    Eigen::MatrixXi F, SF;
    int subdiv = 2;
    igl::readOBJ("../data/cube.obj", V, F);
    //subdivide by 2
    for (int i = 0; i < subdiv; ++i)
    {
        igl::upsample(V, F, SV, SF);
        V = SV;
        F = SF;
    }
    Eigen::Vector3i cluster_size;
    cluster_size << 2, 2, 2;
    add_object_VF(geometry, SV, SF, false, cluster_size);

    Eigen::Vector3d floor_normal;
    Eigen::Vector3d floor_pos;
    Eigen::SparseMatrixd N;
    floor_normal << 0.0, 1.0, 0.0;
    floor_pos << 0.0, -1.0, 0.0;
    //add_plane(floor_normal, floor_pos, geometry);
    add_floor(floor_normal, floor_pos, geometry);

    Visualize::viewer().callback_key_down = key_down_callback;
}
#endif