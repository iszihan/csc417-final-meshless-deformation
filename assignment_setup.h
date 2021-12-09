#ifndef ASSIGNMENT_SETUP_H
#define ASSIGNMENT_SETUP_H

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

//collision detection stuff
#include <collision_detection.h>

//Variable for geometry
Eigen::MatrixXd V, V_skin; //vertices of simulation mesh //this will hold all individual pieces of cloth, I'll load some offsets
Eigen::MatrixXi F, F_skin; //faces of simulation mesh
Eigen::MatrixXd V_sphere, V_sphere_skin; //vertices of simulation mesh //this will hold all individual pieces of cloth, I'll load some offsets
Eigen::MatrixXi F_sphere, F_sphere_skin; //faces of simulation mesh
Eigen::SparseMatrixd N;
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> still_geometry;
std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> moving_geometry;

//material parameters
double mass = 0.1;

//BC
std::vector<unsigned int> fixed_point_indices;
Eigen::SparseMatrixd P;
Eigen::VectorXd x0;

Eigen::SparseMatrixd M; //mass matrix
Eigen::Vector3d center_of_mass; //original center of mass

//scratch memory for assembly
Eigen::VectorXd tmp_qdot;
Eigen::VectorXd tmp_force;
Eigen::VectorXd gravity;
Eigen::VectorXd qtmp;

std::vector<std::pair<Eigen::Vector3d, unsigned int>> spring_points; 
std::vector<std::pair<Eigen::Vector3d, unsigned int>> collision_points;

//collision detection stuff
bool collision_detection_on = false;
bool simulation_pause = true;

//selection spring
double lspring = 0.1;
double k_selected = 1e5;
double k_collision = 1e5;
inline void simulate(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double t)
{
    //std::cout<<"inside simulate"<<std::endl;

    //Interaction spring
    if (! simulation_pause){
        spring_points.clear();
        Eigen::Vector3d mouse;
        Eigen::Vector6d dV_mouse;
        double k_selected_now = (Visualize::is_mouse_dragging() ? k_selected : 0.);
        for (unsigned int pickedi = 0; pickedi < Visualize::picked_vertices().size(); pickedi++)
        {   
            // std::cout<<pickedi<<std::endl;
            Eigen::Vector3d p1 = (q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6);
            Eigen::Vector3d p2 = (q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]);
            // std::cout<<p1<<std::endl;
            // std::cout<<p2<<std::endl;
            //spring_points.push_back(std::make_pair((P.transpose() * q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi]));
            spring_points.push_back(std::make_pair((q + x0).segment<3>(3 * Visualize::picked_vertices()[pickedi]) + Visualize::mouse_drag_world() + Eigen::Vector3d::Constant(1e-6), 3 * Visualize::picked_vertices()[pickedi]));
            
            //TODO: add a dragging handle visualization
            //Visualize::scale_x(1, 2.0);
        }

        //add collision spring points with still geometry
        collision_points.clear();
        Eigen::Vector6d dV_collide;
        for(unsigned int mi=0; mi<moving_geometry.size(); ++mi) {
            for(unsigned int si=0; si<still_geometry.size(); ++si){
                std::cout<<"checking collision"<<std::endl;
                collision_detection(collision_points, mi, si, 
                                    q, still_geometry[si].first, still_geometry[si].second);
            }
        }

        // std::cout<<collision_points.size()<<std::endl;
        // if(collision_points.size()>0){
        //     std::cout<<"colliding points:"<<std::endl;
        //     std::cout<<"p1:"<<std::endl;
        //     std::cout<<moving_geometry[0].first.row(collision_points[0].second)<<std::endl;
        //     std::cout<<"ps:"<<std::endl;
        //     std::cout<<collision_points[0].first<<std::endl;
        //     std::exit(1);
        // }

        //add collision spring points between moving geometries 
        auto force = [&](Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q2, Eigen::Ref<const Eigen::VectorXd> qdot2)
        {
            f = -gravity;

            //dragging force
            // for (unsigned int pickedi = 0; pickedi < spring_points.size(); pickedi++)
            // {
            //     //dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose() * q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            //     dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
            //     std::cout<<"q1:"<<spring_points[pickedi].first<<std::endl;
            //     std::cout<<"q2:"<<(q2 + x0).segment<3>(spring_points[pickedi].second)<<std::endl;
            //     f.segment<3>(3 * Visualize::picked_vertices()[pickedi]) -= dV_mouse.segment<3>(3);
            //     std::cout<<"force:"<<std::endl;
            //     std::cout<<dV_mouse.segment<3>(3)<<std::endl;
            // }

            //collision force
            for (unsigned int ci = 0; ci < collision_points.size(); ci++)
            {
                //dV_spring_particle_particle_dq(dV_mouse, spring_points[pickedi].first, (P.transpose() * q2 + x0).segment<3>(spring_points[pickedi].second), 0.0, k_selected_now);
                dV_spring_particle_particle_dq(dV_collide, collision_points[ci].first, (q2 + x0).segment<3>(3 * collision_points[ci].second), 0.0, k_collision);
                std::cout<<"q1:"<<collision_points[ci].first<<std::endl;
                std::cout<<"q2:"<<(q2 + x0).segment<3>(3 * collision_points[ci].second)<<std::endl;
                f.segment<3>(3 * collision_points[ci].second) += dV_collide.segment<3>(3);
                std::cout<<"added force:"<<std::endl;
                std::cout<<dV_collide.segment<3>(3)<<std::endl;
                std::cout<<"current force:"<<std::endl;
                std::cout<<f.segment<3>(3 * collision_points[ci].second)<<std::endl;
            }
            //std::cout<<"here"<<std::endl;
            //std::cout<<gravity.rows()<<std::endl;
            //std::cout<<P.rows()<<std::endl;
            f = P * f;
            //std::cout<<"here2"<<std::endl;
        };

        meshless_implicit_euler(q, qdot, dt, mass, V, center_of_mass, force, tmp_force);
    }
}

inline void draw(Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, double t)
{
    //update vertex positions using simulation
    Visualize::update_vertex_positions(0, P.transpose() * q + x0);
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
    return false;
}

inline void assignment_setup(int argc, char **argv, Eigen::VectorXd &q, Eigen::VectorXd &qdot)
{
    //load moving geometry data
    igl::readOBJ("../data/cube.obj", V, F);
    //setup simulation
    init_state(q, qdot, V);

    //add geometry to scene
    V_skin = V;
    F_skin = F;
    N.resize(V.rows(), V.rows());
    N.setIdentity();
    Visualize::add_object_to_scene(V, F, V_skin, F_skin, N, Eigen::RowVector3d(244, 165, 130) / 255.);    
    moving_geometry.push_back(std::make_pair(V, F));
    
    //mass matrix
    mass_matrix_particles(M, q, mass);
    if (M.rows() == 0)
    {
        std::cout<<"Mass matrix not implemented, exiting.\n";
        exit(1);
    }

    //constant gravity vector
    gravity.resize(q.rows(), 1);
    dV_cloth_gravity_dq(gravity, M, Eigen::Vector3d(0, -9.8, 0));
    //center of mass 
    center_of_mass = V.colwise().mean();
    // std::cout<<V<<std::endl;
    // std::cout<<center_of_mass<<std::endl;

    // //fix to the floor
    //find_min_vertices(fixed_point_indices, V, 0.001);
    // P.resize(q.rows(), q.rows());
    // P.setIdentity();
    //fixed_point_constraints(P, q.rows(), fixed_point_indices);
    //x0 = q - P.transpose() * P * q; //vector x0 contains position of all fixed nodes, zero for everything else
    // //correct M, q and qdot so they are the right size
    // q = P * q;
    // qdot = P * qdot;
    // M = P * M * P.transpose();

    //not fixed to the floor 
    P.resize(q.rows(), q.rows());
    P.setIdentity();
    x0.resize(q.size());
    x0.setZero();

    //add floor
    Eigen::Vector3d floor_normal;
    Eigen::Vector3d floor_pos;
    floor_normal << 0.7, 0.7, 0.;
    floor_pos << 0., -3.0, 0.;    
    Eigen::MatrixXd V_floor;
    Eigen::MatrixXi F_floor;
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
    for(unsigned int iv=0; iv<V_floor.rows(); ++iv) {
        Eigen::Vector3d rotated = floor_R*V_floor.row(iv).transpose();
        V_floor.row(iv) = (rotated + floor_pos + Eigen::Vector3d(0., -0.01, 0.)).transpose();
    }    
    Visualize::add_object_to_scene(V_floor, F_floor, V_floor, F_floor, N, Eigen::RowVector3d(64,165,130)/255.);
    still_geometry.push_back(std::make_pair(V_floor, F_floor));

    Visualize::viewer().callback_key_down = key_down_callback;
    std::cout<<"finished set up"<<std::endl;
}

#endif
