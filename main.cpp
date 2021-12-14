#include <iostream>
#include <thread>
#include <mutex>
#include <assignment_setup.h>
#include <visualization.h>

//Simulation State
std::mutex mtx;
std::vector<scene_object> geometry;
Eigen::VectorXd q;
Eigen::VectorXd qdot;

double t = 0; //simulation time 
double dt = 0.0001; //time step

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

    //setup
    //assignment_setup(argc, argv, geometry);
    clustering_setup(argc, argv, geometry);

    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    //setup libigl viewer and activate 
    Visualize::setup(q, qdot, true);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();
    return 1; 
}
