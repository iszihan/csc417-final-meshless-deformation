#include <iostream>
#include <thread>
#include <mutex>
#include <visualization.h>
#include <assignment_setup.h>

//Simulation State
std::mutex mtx;
//typedef std::tuple<int, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::SparseMatrixd,
//Eigen::SparseMatrixd, Eigen::Vector3d, Eigen::VectorXd, Eigen::VectorXd, Eigen::SparseMatrixd, Eigen::VectorXd, double, Eigen::Vector3d> scene_object;

typedef std::tuple<int,                           //0 moving or still
	               Eigen::MatrixXd,               //1 V
	               Eigen::MatrixXi,               //2 F
	               Eigen::MatrixXd,               //3 V_skin
	               Eigen::MatrixXi,               //4 F_skin
	               Eigen::SparseMatrixd,          //5 N skinning matrix
	               Eigen::SparseMatrixd,          //6 M -- this is not actually used anywhere?
	               std::vector<Eigen::Vector3d>,  //7 center of mass for clusters
	               Eigen::VectorXd,               //8 q
	               Eigen::VectorXd,               //9 qdot
	               Eigen::SparseMatrixd,          //10 P
	               Eigen::VectorXd,               //11 x0
	               Eigen::VectorXd,               //12 gravity
	               std::vector<std::vector<int>>, //13 clusters of pair of vertex indices and positions
	               std::vector<Eigen::MatrixXd>,  //14 Q = V-center_of_mass for clusters
	               double,                        //15 distance to com
	               Eigen::Vector3d,                //16 com that moves with time
					std::vector<std::vector<int>>   //17  vertex face list
> scene_object;

std::vector<scene_object> geometry;

Eigen::VectorXd q;
Eigen::VectorXd qdot;
//simulation time and time step
double t = 0; //simulation time 
double dt = 0.0001; //time step

//simulation loop
bool simulating = true;
bool simulation_callback() {
    while(simulating){
	 	//simulate_clustering(geometry, dt, t);
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
    assignment_setup(argc, argv, geometry);
    //clustering_setup(argc, argv, geometry);
    
    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();

    //setup libigl viewer and activate 
    Visualize::setup(q, qdot, true);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();
    return 1; 
}
