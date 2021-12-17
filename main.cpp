#include <iostream>
#include <thread>
#include <mutex>
#include <visualization.h>
#include <assignment_setup.h>

//Simulation State
std::mutex mtx;
std::vector<scene_object> geometry;

Eigen::VectorXd q;
Eigen::VectorXd qdot;
//simulation time and time step
double t = 0; //simulation time 
double dt = 0.003; //time step

//simulation loop
bool simulating = true;
bool simulation_callback(Eigen::Vector3i force_setup) {
    
    while(simulating){
	 	simulate(geometry, dt, t, mtx, force_setup);
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
    Eigen::Vector3i force_setup;
    setup(argc, argv, geometry,force_setup, t);
    
    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback, force_setup);
    simulation_thread.detach();

    //setup libigl viewer and activate 
    Visualize::setup(q, qdot, true);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();
    return 1; 
}
