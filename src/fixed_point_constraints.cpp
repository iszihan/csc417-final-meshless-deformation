#include <fixed_point_constraints.h>
#include <algorithm>
#include <set>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {

    // q_size = 3 * n
    // P is a selection matrix of size 3(n-l) x 3n, l is the number of fixed particles ;
    int n = q_size / 3;
    int q_unfixed_size = ( n - indices.size() ) * 3;

    //find the unfixed indices 
    std::vector<unsigned int> indices_all(n);    
    std::vector<int> indices_unfixed;
    std::set<unsigned int> indices_fixed(indices.begin(), indices.end());
    for(int i=0; i<n; ++i){
        indices_all.push_back(i);
        if(indices_fixed.find(i) == indices_fixed.end()){
            indices_unfixed.push_back(i);
        }
    }
    
    // std::cout<<"total # of v :"<<n<<std::endl;
    // std::cout<<"fixed # of v :"<<indices.size()<<std::endl;
    // std::cout<<"unfixed # of v :"<<indices_unfixed.size()<<std::endl;
    P.resize(q_unfixed_size, q_size);
    std::vector<Eigen::Triplet<double>> TripletList ;
    TripletList.reserve(q_unfixed_size);
    for(int i=0; i<indices_unfixed.size(); ++i){
        TripletList.push_back(Eigen::Triplet<double>(i*3,indices_unfixed[i]*3,1));
        TripletList.push_back(Eigen::Triplet<double>(i*3+1,indices_unfixed[i]*3+1,1));
        TripletList.push_back(Eigen::Triplet<double>(i*3+2,indices_unfixed[i]*3+2,1));
    }
    P.setFromTriplets(TripletList.begin(), TripletList.end());
    
}