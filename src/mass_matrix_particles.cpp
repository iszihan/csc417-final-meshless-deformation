#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    //same mass for all particles 
    std::vector<Eigen::Triplet<double>> tripletList ; 
    tripletList.reserve(q.size());
    for (int i=0; i<q.size(); ++i){
         tripletList.push_back(Eigen::Triplet<double>(i,i,mass));
    }
    M.resize(q.size(), q.size());
    M.setFromTriplets(tripletList.begin(),tripletList.end());
}