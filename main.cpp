#include <iostream>
#include <thread>
#include <mutex>
#include <assignment_setup.h>
#include <visualization.h>

//Simulation State
std::vector<Eigen::VectorXd> q_list;
std::vector<Eigen::VectorXd> qdot_list;
std::mutex mtx;
//typedef std::tuple<int, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::SparseMatrixd,
//Eigen::SparseMatrixd, Eigen::Vector3d, Eigen::VectorXd, Eigen::VectorXd, Eigen::SparseMatrixd, Eigen::VectorXd, double, Eigen::Vector3d> scene_object;

typedef std::tuple<int, //moving or still
    Eigen::MatrixXd, //V 
    Eigen::MatrixXi, //F
    Eigen::MatrixXd, //V_skin
    Eigen::MatrixXi, //F_skin
    Eigen::SparseMatrixd, //N skinning matrix
    Eigen::SparseMatrixd, //M
    Eigen::Vector3d, //center of mass
    Eigen::VectorXd, //q
    Eigen::VectorXd, //qdot
    Eigen::SparseMatrixd, //P
    Eigen::VectorXd, //x0
    Eigen::VectorXd, //gravity
    std::vector<std::vector<int>>, //clusters
    Eigen::MatrixXd> scene_object;  //Q = V-center_of_mass

std::vector<scene_object> geometry;

// 0 object type, 0 == plane, 1 == movable object
// 1 V,
// 2 F,
// 3 skinning V,
// 4 skinning F,
// 5 N,
// 6 M,
// 7 com,
// 8 q,
// 9 qdot
// 10 P
// 11 gravity
// 12 radius of com
// 13 current center of mass

Eigen::VectorXd q;
Eigen::VectorXd qdot;
//simulation time and time step
double t = 0; //simulation time 
double dt = 0.00001; //time step

//simulation loop
bool simulating = true;
bool simulation_callback() {
    
    //simulate(q, qdot, dt, t);
    //simulate(q, qdot, dt, t);
    //simulate(q, qdot, dt, t);

    while(simulating) {
		simulate(geometry, dt, t, mtx);
    	t += dt;
    }
    return false;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer) {
    mtx.lock();
    draw(geometry, t);
    mtx.unlock();
    return false;
}

void f(int &a, int b, int c) {
    a = b + c;
}

void g(Eigen::Vector3d &a, Eigen::Vector3d b, Eigen::Vector3d c) {
    a = b + c;
}

template<typename Ret, typename B, typename C>
void h(Ret &&a, B b, C c, void (*func)(Ret, B, C)) {
    func(a,b,c);
}

int main(int argc, char **argv) {
    std::cout<<"Start Meshless Deformation...\n";

    //assignment specific setup
    assignment_setup(argc, argv, q_list, qdot_list, geometry);

    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    //setup libigl viewer and activate 
    Visualize::setup(q, qdot, true);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();
    return 1; 
}
