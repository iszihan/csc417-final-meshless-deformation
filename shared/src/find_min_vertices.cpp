#include <find_min_vertices.h>
#include <iostream>

void find_min_vertices(std::vector<unsigned int> &indices, Eigen::Ref<const Eigen::MatrixXd> V, double tol) {
    
    double min_vertex = V(0,1); 
    
    for(unsigned int vi=0; vi<V.rows(); ++vi) {
        min_vertex = (V(vi,1) < min_vertex ? V(vi,1) : min_vertex);
    }
    std::cout<<min_vertex<<std::endl;

    for(unsigned int vi=0; vi<V.rows(); ++vi) {

        if(std::abs(V(vi,1)-min_vertex) <= tol) {
            indices.push_back(vi);
        }
    }
    // int min_idx = 0;
    // double min_x = V(min_idx,0); 
    // double min_y = V(min_idx,1); 
    // double min_z = V(min_idx,2); 
    // for(unsigned int vi=0; vi<V.rows(); ++vi) {
    //     if(V.row(vi)(0) < min_x && V.row(vi)(1) < min_y && V.row(vi)(2) < min_z){
    //         min_idx = vi;
    //         min_x = V(min_idx,0); 
    //         min_y = V(min_idx,1); 
    //         min_z = V(min_idx,2); 
    //     }
    // }
    // indices.push_back(min_idx);
}