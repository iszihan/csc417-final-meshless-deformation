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

//data structures for geometry
// typedef std::tuple<int,                           //0 moving or still
//     Eigen::MatrixXd,               //1 V
//     Eigen::MatrixXi,               //2 F
//     Eigen::MatrixXd,               //3 V_skin
//     Eigen::MatrixXi,               //4 F_skin
//     Eigen::SparseMatrixd,          //5 N skinning matrix
//     Eigen::SparseMatrixd,          //6 M -- this is not actually used anywhere?
//     std::vector<Eigen::Vector3d>,  //7 center of mass for clusters
//     Eigen::VectorXd,               //8 q
//     Eigen::VectorXd,               //9 qdot
//     Eigen::SparseMatrixd,          //10 P
//     Eigen::VectorXd,               //11 x0
//     Eigen::VectorXd,               //12 gravity
//     std::vector<std::vector<int>>, //13 clusters of pair of vertex indices and positions
//     std::vector<Eigen::MatrixXd>,  //14 Q = V-center_of_mass for clusters
//     double,                        //15 distance to com
//     Eigen::Vector3d,                //16 com that moves with time
//     std::vector<std::vector<int>>,   //17  vertex face list
//     std::unordered_map <Eigen::Vector3d, Bounding_box, Spatial_hash_fn>, // 18 spacial hash table
//     std::vector<int>				// 19 occupied list
// > scene_object;

std::string data_paths[3] = {"../data/cube.obj",
                             "../data/coarse_bunny2.obj",
                             "../data/sphere.obj"};

//material parameters
double mass = 1.0;

//scratch memory for assembly -- do we need these beside tmp_force?
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
std::vector<std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>>> collision_points_list; // 0=point on the surface of collision target, 1=vertex normal, 2=vertex index

//integration method
// 0 - rigid
// 1 - linear 
// 2 - quadratic
int method = 0;
int which_setup = 0;

//collision detection stuff
bool collision_detection_on = false;
bool simulation_pause = true;
bool drop_key_pressed = false;
int item_type = 0;

//selection spring
double lspring = 0.1;
double k_selected = 1e7;
double k_collision = 1e4;
inline void add_object(std::vector<scene_object>& geometry, std::string file_path, 
                       Eigen::Vector3d position, Eigen::Vector3d scale, int subdiv,
                       Eigen::Vector3i clusters, double fixed_tol, bool fixed)
{
    Eigen::VectorXd q;
    Eigen::VectorXd qdot;
    Eigen::VectorXd gravity;
    Eigen::SparseMatrixd M;
    Eigen::MatrixXd V, V_skin; //vertices of simulation mesh
    Eigen::MatrixXi F, F_skin; //faces of simulation mesh
    Eigen::SparseMatrixd N;
    Eigen::SparseMatrixd P;
    Eigen::VectorXd x0;
    Eigen::MatrixXd SV;
    Eigen::MatrixXi SF;

    //load in file
    igl::readOBJ(file_path, V, F);
    for (int i = 0; i < subdiv; ++i)
    {
        igl::upsample(V, F, SV, SF);
        V = SV;
        F = SF;
    }
    //translate to origin 
    Eigen::Vector3d com = V.colwise().mean();
    V = V.rowwise() - com.transpose();
    //scale
    V.array().rowwise() *= scale.transpose().array();
    //translate to target position
    V = V.rowwise() + position.transpose();
    init_state(q, qdot, V);

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
        find_min_vertices(fixed_point_indices, V, fixed_tol);
        P.resize(q.rows(), q.rows());
        P.setIdentity();
        fixed_point_constraints(P, q.rows(), fixed_point_indices);
        x0 = q - P.transpose() * P * q; //vector x0 contains position of all fixed nodes, zero for everything else
    }
    else
    {
        //not fixed to the floor
        P.resize(q.rows(), q.rows());
        P.setIdentity();
        x0.resize(q.size());
        x0.setZero();
    }

    //form clusters
    std::vector<Eigen::Vector3d> centers_of_mass;
    std::vector<Eigen::MatrixXd> Qs;
    std::vector<Eigen::Matrix3d> Sps; //for plasticity
    std::vector<std::vector<int>> vertex_clusters;
    Eigen::VectorXi vertex_cluster_counts(V.rows()); //how many cluster each vertex belongs to 
    //add the global cluster
    std::vector<int> all_vertices;
    for(int i=0; i<V.rows(); ++i){
        all_vertices.push_back(i);
        vertex_cluster_counts(i) = 1;
    }
    Eigen::Vector3d center_of_mass = V.colwise().mean();
    Eigen::MatrixXd Q = V.rowwise() - center_of_mass.transpose();
    Sps.push_back(Eigen::Matrix3d::Identity());
    centers_of_mass.push_back(center_of_mass);
    vertex_clusters.push_back(all_vertices);
    Qs.push_back(Q);

    //add sub-clusters if any
    if(clusters(0)>1 || clusters(1)>1 || clusters(2)>1){
        double overlap = 1.2;
        Eigen::MatrixXd V_cluster;
        Eigen::Vector3d V_max = V.colwise().maxCoeff();
        Eigen::Vector3d V_min = V.colwise().minCoeff();
        double dx = (V_max(0) - V_min(0)) / clusters(0);
        double dy = (V_max(1) - V_min(1)) / clusters(1);
        double dz = (V_max(2) - V_min(2)) / clusters(2);
        for (int ix = 0; ix < clusters(0); ++ix)
        {
            double x_min = V_min(0) + dx * (ix - overlap);
            double x_max = V_min(0) + dx * (ix + overlap);
            for (int iy = 0; iy < clusters(1); ++iy)
            {
                double y_min = V_min(1) + dy * (iy - overlap);
                double y_max = V_min(1) + dy * (iy + overlap);
                for (int iz = 0; iz < clusters(2); ++iz)
                {
                    double z_min = V_min(2) + dz * (iz - overlap);
                    double z_max = V_min(2) + dz * (iz + overlap);

                    std::vector<int> curr_cluster;
                    curr_cluster.clear();
                    V_cluster.setZero();
                    double tol = 1e-3;
                    for (int iv = 0; iv < V.rows(); ++iv)
                    {
                        //put vertices to current cluster
                        if (V.row(iv)(0) >= x_min - tol && V.row(iv)(0) <= x_max + tol && V.row(iv)(1) >= y_min - tol && V.row(iv)(1) <= y_max + tol && V.row(iv)(2) >= z_min - tol && V.row(iv)(2) <= z_max + tol)
                        {
                            curr_cluster.push_back(iv);
                            V_cluster.conservativeResize(curr_cluster.size(),3);
                            V_cluster.row(curr_cluster.size()-1) = V.row(iv);
                            vertex_cluster_counts(iv) += 1;
                        }
                    }
                    vertex_clusters.push_back(curr_cluster);
                    Eigen::Vector3d curr_com = V_cluster.colwise().mean();
                    Eigen::MatrixXd curr_Q = V_cluster.rowwise() - curr_com.transpose();
                    centers_of_mass.push_back(curr_com);
                    Qs.push_back(curr_Q);
                    Sps.push_back(Eigen::Matrix3d::Identity());
                    // std::cout<<"curr cluster size:"<<curr_cluster.size()<<std::endl;
                    //std::cout<<"vertex cluster count:"<<vertex_cluster_counts<<std::endl;
                }
            }
        }
    }

    //append everything to the corresponding list
    scene_object one_geometry;
    std::get<0>(one_geometry) = 1;
    std::get<1>(one_geometry) = V;
    std::get<2>(one_geometry) = F;
    std::get<3>(one_geometry) = V;
    std::get<4>(one_geometry) = F;
    std::get<5>(one_geometry) = N;
    std::get<6>(one_geometry) = M;
    std::get<7>(one_geometry) = centers_of_mass;
    std::get<8>(one_geometry) = q;
    std::get<9>(one_geometry) = qdot;
    std::get<10>(one_geometry) = P;
    std::get<11>(one_geometry) = x0;
    std::get<12>(one_geometry) = gravity;
    std::get<13>(one_geometry) = vertex_clusters;
    std::get<14>(one_geometry) = Qs;
    double radius = ((V.rowwise() - center_of_mass.transpose()).rowwise().norm()).maxCoeff() * 1.5;
    std::get<15>(one_geometry) = radius;
    std::get<16>(one_geometry) = center_of_mass;
    std::vector<std::vector<int>> v2f;
    compute_vertex_face_list(V, F, v2f);
    std::get<17>(one_geometry) = v2f;
    std::get<20>(one_geometry) = vertex_cluster_counts;
    std::get<21>(one_geometry) = Sps;
    geometry.push_back(one_geometry);
    Visualize::add_object_to_scene(V, F, V, F, N, Eigen::RowVector3d(244, 165, 130) / 255.);
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
    float angle = std::acos(n0.dot(floor_normal));
    Eigen::Vector3d axis = n0.cross(floor_normal);
    Eigen::Matrix3d floor_R = Eigen::AngleAxisd(angle, axis).matrix();
    //translate plane       
    for (unsigned int iv = 0; iv < V_floor.rows(); ++iv)
    {
        Eigen::Vector3d rotated = floor_R * V_floor.row(iv).transpose();
        V_floor.row(iv) = (rotated + floor_pos + Eigen::Vector3d(0., -0.01, 0.)).transpose();
    }

    Visualize::add_object_to_scene(V_floor, F_floor, V_floor, F_floor, N, Eigen::RowVector3d(64, 165, 130) / 255.);

    Eigen::VectorXd q(V_floor.rows()*3);
    q.setZero();
    scene_object one_geometry;
    Eigen::SparseMatrixd M;
    std::get<0>(one_geometry) = 0;
    std::get<1>(one_geometry) = V_floor;
    std::get<2>(one_geometry) = F_floor;
    std::get<3>(one_geometry) = V_floor;
    std::get<4>(one_geometry) = F_floor;
    std::get<5>(one_geometry) = N;
    std::get<6>(one_geometry) = M;
    std::get<7>(one_geometry).push_back(Eigen::Vector3d::Zero());
    std::get<8>(one_geometry) = q;
    std::get<9>(one_geometry) = q;
    geometry.push_back(one_geometry);
}
inline void add_floor_alt(Eigen::Vector3d floor_normal, Eigen::Vector3d floor_pos, std::vector<scene_object>& geometry)
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
    float angle = std::acos(n0.dot(floor_normal));
    Eigen::Vector3d axis = n0.cross(floor_normal);
    Eigen::Matrix3d floor_R = Eigen::AngleAxisd(angle, axis).matrix();
    //translate plane       
    for (unsigned int iv = 0; iv < V_floor.rows(); ++iv)
    {
        Eigen::Vector3d rotated = floor_R * V_floor.row(iv).transpose();
        V_floor.row(iv) = (rotated + floor_pos + Eigen::Vector3d(0., -0.01, 0.)).transpose();
    }

    Visualize::add_object_to_scene(V_floor, F_floor, V_floor, F_floor, N, Eigen::RowVector3d(64, 165, 130) / 255.);

    Eigen::VectorXd q(V_floor.rows() * 3);
    Eigen::VectorXd qdot(V_floor.rows() * 3);
    q.setZero();
    init_state(q, qdot, V_floor);
    scene_object one_geometry;
    Eigen::SparseMatrixd M;
    Eigen::Vector3d curr_com = V_floor.colwise().mean();
	
    std::get<0>(one_geometry) = -1;
    std::get<1>(one_geometry) = V_floor;
    std::get<2>(one_geometry) = F_floor;
    std::get<3>(one_geometry) = V_floor;
    std::get<4>(one_geometry) = F_floor;
    std::get<5>(one_geometry) = N;
    std::get<6>(one_geometry) = M;
    std::get<7>(one_geometry).push_back(Eigen::Vector3d::Zero());
    std::get<8>(one_geometry) = q;
    std::get<9>(one_geometry) = q;
    std::get<16>(one_geometry) = curr_com;
    geometry.push_back(one_geometry);
}
inline void simulate(std::vector<scene_object> &geometry, double dt, double t, std::mutex &mtx, Eigen::Vector4i force_setup)
{
    //std::cout<<"inside simulate"<<std::endl;
    bool add_gravity = (force_setup(0) == 1);
    bool add_dragging = (force_setup(1) == 1);
    bool add_collision = (force_setup(2) == 1);
    bool add_collision_optimization = (force_setup(3) == 1);

    //Interaction spring
    if (!simulation_pause)
    {

        spring_points_list.clear();
        collision_points_list.clear();
        for (int i = 0; i < geometry.size(); i++)
        {
            std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
            std::vector < std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>> collision_points_tmp;
            spring_points_list.push_back(spring_points_tmp);
            collision_points_list.push_back(collision_points_tmp);
            if (add_collision_optimization) {
                construct_spatial_hash_table(geometry.at(i));
            }
        }

        if(add_collision){
            //add collision spring points with other geometry
            Eigen::Vector6d dV_collide;
            for (int mi = 0; mi < geometry.size(); ++mi)
            {
                std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points_tmp;
                std::vector < std::tuple<Eigen::Vector3d, Eigen::Vector3d, unsigned int, unsigned int, unsigned int>> collision_points_tmp;
                scene_object moving_object = geometry.at(mi);
                if (std::get<0>(moving_object) >= 1) //only check if it's static
                {
                    for (unsigned int si = 0; si < geometry.size(); ++si)
                    {   
                        if (si != mi) //do not check against itself
                        {
                            //std::cout<<si<<std::endl;
                            scene_object collision_target = geometry.at(si); 
                            // only do collision compute if they overlaps
                            if (std::get<0>(collision_target) >= 1) {
                                if (precomputation(moving_object, collision_target)) {
                                    // only have collision detection with planes for now others can be implemented later
                                    if (add_collision_optimization) {
                                        collision_detection_with_optimization(collision_points_list.at(mi), mi, si, moving_object, collision_target);
                                    }
                                    else
                                    {
                                        collision_detection(collision_points_list.at(mi), mi, si, moving_object, collision_target);
                                    }
                                    //std::cout << float(clock() - begin_time) / CLOCKS_PER_SEC << std::endl;
                                }
                            } else
                            {
                                if (add_collision_optimization) {
                                    collision_detection_with_optimization(collision_points_list.at(mi), mi, si, moving_object, collision_target);
                                }
                                else
                                {
                                    collision_detection(collision_points_list.at(mi), mi, si, moving_object, collision_target);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        Eigen::Vector6d dV_mouse;
        double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
        if(add_dragging){
            for (int object_id = 0; object_id < geometry.size(); object_id++)
            {
                if (std::get<0>(geometry.at(object_id)) > 0)  // only consider the movable objects to save compute
                {
                    std::cout<<Visualize::picked_vertices().size()<<std::endl;
                    for (unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++)
                    {   
                        std::cout<<Visualize::picked_vertices()[pickedi].first<<std::endl;
                        if(object_id == Visualize::picked_vertices()[pickedi].first){
                            Eigen::VectorXd q = std::get<8>(geometry.at(object_id));
                            Eigen::Vector3d p1 = q.segment<3>(3 * Visualize::picked_vertices()[pickedi].second) + Visualize::mouse_drag_world() * 20.f + Eigen::Vector3d::Constant(1e-6);
                            Eigen::Vector3d p2 = q.segment<3>(3 * Visualize::picked_vertices()[pickedi].second);
                            spring_points_list.at(object_id).push_back(std::make_pair(q.segment<3>(3 * Visualize::picked_vertices()[pickedi].second) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi].second));
                        }
                    }
                }
            }
        }

        for (int i = 0; i < geometry.size(); i++)
        {
            scene_object current_object = geometry.at(i);
            if (std::get<0>(current_object) > 0)
            {
                auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2)
                {
                    Eigen::SparseMatrixd P = std::get<10>(current_object);
                    std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points = spring_points_list.at(i);
                    std::cout<<"spring points size:"<<spring_points.size()<<std::endl;
                    
                    //gravity
                    f = -std::get<12>(current_object);
                    if (!add_gravity){
                        f.setZero();
                    }

                    //dragging force
                    if(add_dragging){
                        for (unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++)
                        {
                           dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, q2.segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                           f.segment<3>(3 * Visualize::picked_vertices()[pickedi].second) -= dV_mouse.segment<3>(3);
                        }
                    }

                    //collision force
                    if(add_collision){
                        // if(collision_points_list.at(i).size()>0){
                        //     std::cout<<"collision size here:"<<collision_points_list.at(i).size()<<std::endl;
                        // }                         
                        for (unsigned int ci = 0; ci < collision_points_list.at(i).size(); ci++)
                        {
                            //std::cout<<"adding collision "<<ci<<std::endl;
                            Eigen::Vector3d repulsive_force, pt, pt_projected, target_dir, ptdot, pt2dot;
                            //current colliding vertex projected onto the colliding object
                            pt_projected = std::get<0>(collision_points_list.at(i).at(ci));
                            //current colliding vertex
                            pt = std::get<8>(current_object).segment<3>(std::get<2>(collision_points_list.at(i).at(ci)) * 3);
                            //current colliding vertex velocity
                            ptdot = std::get<9>(current_object).segment<3>(std::get<2>(collision_points_list.at(i).at(ci)) * 3);
                            //the other colliding vertex id
                            int obj2_vertex_id_velocity = std::get<3>(collision_points_list.at(i).at(ci));
                            //the other colliding vertex's velocity
                            pt2dot = std::get<9>(geometry.at(std::get<4>(collision_points_list.at(i).at(ci)))).segment<3>(obj2_vertex_id_velocity * 3);
                            target_dir = std::get<1>(collision_points_list.at(i).at(ci)).normalized(); // this is made negative, so that it becomes the direction the repulsive force should go to
                            double d = abs((pt - pt_projected).dot(-target_dir));
                            double force_magnitude = d * k_collision - 5000 * (ptdot.dot(-target_dir) - pt2dot.dot(-target_dir));
                            force_magnitude = d * k_collision;
                            repulsive_force = force_magnitude * (-target_dir); // this will be in the direction where obj1 will be bounced to
                            f.segment<3>(3 * std::get<2>(collision_points_list.at(i).at(ci))) += repulsive_force;
                            
                            // std::cout<<"the other object type: "<<std::get<9>(geometry.at(std::get<4>(collision_points_list.at(i).at(ci))))<<std::endl;
                            // std::cout<<std::get<2>(collision_points_list.at(i).at(ci))<<std::endl;
                            // std::cout<<"ptdot:"<<ptdot<<std::endl;
                            // std::cout<<"pt2dot:"<<pt2dot<<std::endl;
                            // std::cout<<"pt projected:"<<pt_projected<<std::endl;
                            // std::cout<<"pt:"<<pt<<std::endl;
                            // std::cout<<"force mag:"<<force_magnitude<<std::endl;
                            // std::cout<<"k_collision"<<k_collision<<std::endl;
                            // std::cout<<"d * k_collision"<<d * k_collision<<std::endl;
                            // std::cout<<"d"<<d<<std::endl;
                            // std::cout<<repulsive_force<<std::endl;
                        }
                        // if(collision_points_list.at(i).size()>0){
                        //     std::cout<<"collision size here:"<<collision_points_list.at(i).size()<<std::endl;
                        //     std::cout<<f.norm()<<std::endl;
                        // }
                        // std::cout<<"force:"<<f.norm()<<std::endl;
                    }
                };

                Eigen::Vector3d comt;
                q_tmp = std::get<8>(current_object);
                qdot_tmp = std::get<9>(current_object);
                meshless_implicit_euler(q_tmp, qdot_tmp, dt, mass, 
                                        std::get<10>(current_object), std::get<11>(current_object), 
                                        std::get<7>(current_object), std::get<14>(current_object), 
                                        std::get<13>(current_object), std::get<20>(current_object),
                                        method, force, tmp_force, comt, std::get<21>(current_object));
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
            add_object(geometry, data_paths[item_placement], pos, Eigen::Vector3d(1.0, 1.0, 1.0), 0, Eigen::Vector3i(1, 1, 1), 0.001, false);
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
            Visualize::update_vertex_positions(i, std::get<8>(current_object));
        }
    }
}

bool key_down_callback(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
    if (key == '4')
    {
        std::cout << "toggle deformation to linear \n";
        method = 1;
    }
    if (key == '5')
    {
        std::cout << "toggle deformation to quadratic \n";
        method = 2;
    }
    if (key == '6')
    {
        std::cout << "toggle deformation to rigid only \n";
        method = 0;
    }
    if (key == 'T')
    {
        std::cout << "toggle deformation to plasticity updates \n";
        method = 3;
    }

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

inline void setup(int argc, char **argv, std::vector<scene_object> &geometry, Eigen::Vector4i& force_setup, double& t)
{    

    //things that need to be tweaked
    k_selected = 1e4;
    k_collision = 1e5;
    t = 0.01;

    if(argc>1){
        which_setup = std::stoi(std::string(argv[1]));
    }
    if(which_setup == 0){
        //bunny collision setup
        k_selected = 1e4;
        k_collision = 1e5;
        t = 0.01;
        add_object(geometry, data_paths[1], Eigen::Vector3d(0.0,8.0,-1.0), Eigen::Vector3d(1.0,1.0,1.0), 0, Eigen::Vector3i(1,1,1), 0.001, false);
        add_object(geometry, data_paths[1], Eigen::Vector3d(0.0, 3.0, -1.0), Eigen::Vector3d(1.0, 1.0, 1.0), 0, Eigen::Vector3i(1, 1, 1), 0.001, false);
        add_floor_alt(Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d(0.0,-1.0,0.0), geometry);
        add_floor_alt(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(-17.0, 8, 0.0), geometry);
        add_floor_alt(Eigen::Vector3d(-1.0, 0.0, 0.0), Eigen::Vector3d(17.0, 8, 0.0), geometry);
        add_floor_alt(Eigen::Vector3d(0, 0.0, -1.0), Eigen::Vector3d(0.0, 8, 17.0), geometry);
        add_floor_alt(Eigen::Vector3d(0, 0.0, 1.0), Eigen::Vector3d(0.0, 8, -17.0), geometry);
        force_setup<<1,0,1,1;
    }else if(which_setup == 1){
        //cubes collision setup
        add_object(geometry, data_paths[0], Eigen::Vector3d(1.0, 4.0, 1.0), Eigen::Vector3d(1.0, 1.0, 1.0), 0, Eigen::Vector3i(1, 1, 1), 0.001, false);
        add_object(geometry, data_paths[0], Eigen::Vector3d(1.0, 2.0, 1.0), Eigen::Vector3d(1.0, 1.0, 1.0), 0, Eigen::Vector3i(1, 1, 1), 0.001, false);
        add_object(geometry, data_paths[0], Eigen::Vector3d(0.0,2.0,0.0), Eigen::Vector3d(1.0,1.0,1.0), 0, Eigen::Vector3i(1,1,1), 0.001, false);
        add_floor_alt(Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d(0.0,-1.0,0.0), geometry);
        add_floor_alt(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(-17.0, 8, 0.0), geometry);
        add_floor_alt(Eigen::Vector3d(-1.0, 0.0, 0.0), Eigen::Vector3d(17.0, 8, 0.0), geometry);
        add_floor_alt(Eigen::Vector3d(0, 0.0, -1.0), Eigen::Vector3d(0.0, 8, 17.0), geometry);
        add_floor_alt(Eigen::Vector3d(0, 0.0, 1.0), Eigen::Vector3d(0.0, 8, -17.0), geometry);
        force_setup<<1,0,1,1;
    }else if(which_setup == 2){
        //cube clustering comparison setup
        //1. fixed to the floor dragging scenario to demonstrate clustering
        k_collision = 1e5;
        k_selected = 2e7;
        t = 0.002;
        add_object(geometry, "../data/cube.obj", Eigen::Vector3d(-4.0,3.0,0.0), Eigen::Vector3d(1.0,2.0,1.0), 5, Eigen::Vector3i(1,1,1), 0.001, true);
        add_object(geometry, "../data/cube.obj", Eigen::Vector3d(0.0,3.0,0.0), Eigen::Vector3d(1.0,2.0,1.0), 5, Eigen::Vector3i(1,5,1), 0.001, true);
        add_object(geometry, "../data/cube.obj", Eigen::Vector3d(4.0,3.0,0.0), Eigen::Vector3d(1.0,2.0,1.0), 5, Eigen::Vector3i(1,10,1), 0.001, true);
        add_floor(Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d(0.0,1.0,0.0), geometry);
        force_setup<<0,1,0;
    }else if(which_setup == 3){
        //bunny fixed to the floor setup        
        k_selected = 1e7;
        k_collision = 1e7;
        t = 0.001;
        add_object(geometry, "../data/coarse_bunny2.obj", Eigen::Vector3d(0.0,5.0,0.0), Eigen::Vector3d(1.0,1.0,1.0), 0, Eigen::Vector3i(1,1,1), 0.5, true);
        force_setup<<0,1,0;
        //bunny scenario
        add_floor(Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d(0.0,1.0,0.0), geometry);
    }
    Visualize::viewer().callback_key_down = key_down_callback;
}

#endif