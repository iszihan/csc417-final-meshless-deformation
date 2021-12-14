#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H
#include <tuple>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readOFF.h>
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

//collision detection stuff
#include <collision_detection.h>

//Variable for geometry
typedef std::tuple<int, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::SparseMatrixd,
Eigen::SparseMatrixd, Eigen::Vector3d, Eigen::VectorXd, Eigen::VectorXd, Eigen::SparseMatrixd, Eigen::VectorXd, double, Eigen::Vector3d> scene_object;

// export the obj file with these modifications: http://www.opengl-tutorial.org/beginners-tutorials/tutorial-7-model-loading/

std::string data_paths[3] = { "../data/cube.obj", "../data/coarse_bunny2.obj",
                             "../data/cube.obj"};

std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> still_geometry;
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> moving_geometry;
std::vector<Eigen::VectorXd> force_list;
//material parameters
double mass = 0.1;

//BC
std::vector<unsigned int> fixed_point_indices;
Eigen::VectorXd x0;

//scratch memory for assembly
Eigen::VectorXd qdot_tmp;
Eigen::VectorXd tmp_force;
Eigen::VectorXd q_tmp;
Eigen::MatrixXd V_tmp;
Eigen::MatrixXi F_tmp;
Eigen::Vector3d com_tmp0, com_tmpt;
int item_placement = -1;


std::vector<std::vector<std::pair<Eigen::Vector3d, unsigned int>>> spring_points_list;
std::vector<std::vector<std::pair<Eigen::Vector3d, unsigned int>>> collision_points_list;

//collision detection stuff
bool collision_detection_on = false;
bool simulation_pause = true;
bool drop_key_pressed = false;
int item_type = 0;
//selection spring
double lspring = 0.1;
double k_selected = 1e5;
double k_collision = 1e5;

inline void add_object(std::vector<scene_object> &geometry, std::string file_path, Eigen::Vector3d position)
{
    Eigen::VectorXd q;
    Eigen::VectorXd qdot;
    Eigen::VectorXd gravity;
    Eigen::SparseMatrixd M;
    Eigen::Vector3d center_of_mass;
    Eigen::MatrixXd V, V_skin; //vertices of simulation mesh //this will hold all individual pieces of cloth, I'll load some offsets
    Eigen::MatrixXi F, F_skin; //faces of simulation mesh
    Eigen::SparseMatrixd N;
    Eigen::SparseMatrixd P;
    igl::readOBJ(file_path, V, F);

    init_state(q, qdot, V);
	for (int vi = 0; vi < q.size()/3; vi++)
	{
        q.segment<3>(vi * 3) += position;
	}
	
    center_of_mass = V.colwise().mean().transpose();
    //add geometry to scene
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

    //constant gravity vector
    gravity.resize(q.rows(), 1);
    dV_cloth_gravity_dq(gravity, M, Eigen::Vector3d(0, -9.8, 0));
    //center of mass
    //not fixed to the floor 
    P.resize(q.rows(), q.rows());
    P.setIdentity();
    x0.resize(q.size());
    x0.setZero();

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
    std::get<11>(one_geometry) = gravity;
    std::get<12>(one_geometry) = (((-V).rowwise() + center_of_mass.transpose()).rowwise().norm()).maxCoeff() * 1.5;
    std::get<13>(one_geometry) = center_of_mass;



	

    geometry.push_back(one_geometry);
    Visualize::add_object_to_scene(V, F, V, F, N, Eigen::RowVector3d(244, 165, 130) / 255.);
}
inline void simulate(std::vector<scene_object> &geometry, double dt, double t, std::mutex & mtx)
{
    //std::cout<<"inside simulate"<<std::endl;
    force_list.clear();
    //Interaction spring
    if (!simulation_pause) {
    	
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
        for (int mi = 0; mi < geometry.size(); ++mi) {
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> collision_points_tmp;
        	
            scene_object moving_object = geometry.at(mi);
        	// only check if the object is not static
        	if (std::get<0>(moving_object) >= 1)
        	{
        		// check it against each other objects
                for (unsigned int si = 0; si < geometry.size(); ++si) {
                    // do not check collision against itself
                    if (si != mi) {
                    	// only consider if it is ball-parked to be touching
                        scene_object collision_target = geometry.at(si);
                    	if (precomputation(moving_object, collision_target)){
                    		// only have collision detection with planes for now others can be implemented later
                            if (std::get<0>(collision_target) == 0) {
                                collision_detection(collision_points_list.at(mi), mi, si,
                                    std::get<8>(moving_object), std::get<1>(collision_target), std::get<2>(collision_target));
                            }
                        }
                    }
                }
        	}
        }

        Eigen::Vector3d mouse;
        Eigen::Vector6d dV_mouse;
        double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
        for (int object_id = 0; object_id < geometry.size(); object_id++) {
        	// only consider the movable objects to save compute
            if (std::get<0>(geometry.at(object_id)) > 0) {
                for (unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++)
                {
                    Eigen::VectorXd q = std::get<8>(geometry.at(object_id));
                    // std::cout<<pickedi<<std::endl;
                    Eigen::Vector3d p1 = (q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6);
                    Eigen::Vector3d p2 = (q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]);
                    // std::cout<<p1<<std::endl;
                    // std::cout<<p2<<std::endl;
                    //spring_points.push_back(std::make_pair((P.transpose() * q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi]));
                    spring_points_list.at(object_id).push_back(std::make_pair((q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi]));
                        //TODO: add a dragging handle visualization
                        //Visualize::scale_x(1, 2.0);
                }
            }
        }
    	
    	
        //for (unsigned int mi = 0; mi < moving_geometry.size(); ++mi) {
        //    for (unsigned int si = 0; si < still_geometry.size(); ++si) {
        //        std::cout << "checking collision" << std::endl;
        //        collision_detection(collision_points, mi, si,
        //            q, still_geometry[si].first, still_geometry[si].second);
        //    }
        //}

        // std::cout<<collision_points.size()<<std::endl;
        // if(collision_points.size()>0){
        //     std::cout<<"colliding points:"<<std::endl;
        //     std::cout<<"p1:"<<std::endl;
        //     std::cout<<moving_geometry[0].first.row(collision_points[0].second)<<std::endl;
        //     std::cout<<"ps:"<<std::endl;
        //     std::cout<<collision_points[0].first<<std::endl;
        //     std::exit(1);
        // }
        for (int i = 0; i < geometry.size(); i++) {
        //add collision spring points between moving geometries
            scene_object current_object = geometry.at(i);
            if (std::get<0>(current_object) > 0) {
                auto force = [&](Eigen::VectorXd& f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2)
                {
                    Eigen::SparseMatrixd P = std::get<10>(current_object);
                    std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points = spring_points_list.at(i);
                    f = -std::get<11>(current_object);
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
                        dV_spring_particle_particle_dq(dV_collide, collision_points_list.at(i)[ci].first, (q2 + x0).segment<3>(3 * collision_points_list.at(i)[ci].second), 0.0, k_collision);
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

                q_tmp = std::get<8>(current_object);
                qdot_tmp = std::get<9>(current_object);
                V_tmp = std::get<1>(current_object);
                F_tmp = std::get<2>(current_object);
                com_tmp0 = std::get<7>(current_object);

                
                meshless_implicit_euler(q_tmp, qdot_tmp, dt, mass, V_tmp, com_tmp0, force, tmp_force, com_tmpt);
                std::get<8>(geometry.at(i)) = q_tmp;
                std::get<9>(geometry.at(i)) = qdot_tmp;
                std::get<13>(geometry.at(i)) = com_tmpt;
            }
		}
        if (item_placement >= 0)
        {
            mtx.lock();
            Eigen::Vector3d pos = Visualize::mouse_world();
            //pos(1) = std::min(pos(1), 3.0);
            add_object(geometry, data_paths[item_placement], pos);
            item_placement = -1;
            mtx.unlock();
        }
    }
}

inline void draw(std::vector<scene_object> geometry, double t)
{
    //update vertex positions using simulation
    for (int i = 0; i < geometry.size(); i++) {
        scene_object current_object = geometry.at(i);
        if (std::get<0>(current_object) > 0) {
		   Visualize::update_vertex_positions(i, std::get<10>(current_object).transpose() * std::get<8>(current_object) + x0);
		}
    }
}

bool key_down_callback(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == 'N')
    {
        std::cout << "toggle integrators \n";
        //fully_implicit = !fully_implicit;
    }
    if (key == 'C')
    {
        collision_detection_on = !collision_detection_on;
        Visualize::set_visible(1, collision_detection_on);
    }
    if (key == 'Q'){
        std::exit(1);
    }
    if (key == 'P'){
        simulation_pause = !simulation_pause;
    }

	if ((key == '1' || key == '2' || key == '3') && item_placement == -1)
	{
		if (key == '1'){
            item_placement = 0;
        } else if (key == '2') {
            item_placement = 1;
        } else if (key == '3') {
            item_placement = 2;
        }
        
        std::cout << "\n gonna be dropping item " << item_placement << "\n";

	}
	
    return false;
}
inline void add_plane(Eigen::Vector3d floor_normal, Eigen::Vector3d floor_pos, std::vector<scene_object> & scene_object_list)
{
    Eigen::MatrixXd V_floor;
    Eigen::MatrixXi F_floor;
    Eigen::SparseMatrixd N;
    igl::readOBJ("../data/plane.obj", V_floor, F_floor);
    V_floor *= 10.0; //make it bigger
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
    for (unsigned int iv = 0; iv < V_floor.rows(); ++iv) {
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
    scene_object_list.push_back(one_geometry);
}
inline void assignment_setup(int argc, char **argv, std::vector<Eigen::VectorXd> &q_list, std::vector<Eigen::VectorXd>&qdot_list, std::vector<scene_object> &geometry)
{
	// add just a cube to the scene
    add_object(geometry, data_paths[0], Eigen::Vector3d::Zero());
	
    //add floor
    Eigen::Vector3d floor_normal;
    Eigen::Vector3d floor_pos;
    floor_normal << 0, 0.7, 0.;
    floor_pos << 0., -3.0, 0.;
    add_plane(floor_normal, floor_pos, geometry);
    
    Visualize::viewer().callback_key_down = key_down_callback;
    std::cout<<"finished set up"<<std::endl;
    std::cout << "there is a total of " << geometry.size() << " pieces of geometry \n";
}

#endif
