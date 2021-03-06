#include <pick_nearest_vertices.h>
#include <iostream>
#include <igl/Hit.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/unproject.h>

bool pick_nearest_vertices(std::vector<std::pair<unsigned int, unsigned int>> &verts, Eigen::Ref<const Eigen::Vector3d> win, 
                           Eigen::Ref<const Eigen::Matrix44f> view, Eigen::Ref<const Eigen::Matrix44f> proj, Eigen::Vector4f viewport,
                           Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double radius, int gi) {

    //verts.clear();

    //based on ray shooting
    //source, destination and direction in world
     Eigen::Vector3f start,dir;
     Eigen::Vector3f win_0(win(0), win(1), win(2));
     Eigen::Vector3f win_1(win(0), win(1), 1.);

     igl::unproject(win_0,view,proj,viewport,start);
     igl::unproject(win_1,view,proj,viewport,dir);
     dir -= start; 

     igl::Hit hit;

     const auto & shoot_ray = [&V,&F](const Eigen::Vector3f& s, const Eigen::Vector3f& dir, igl::Hit& hit)->bool
     {
         std::vector<igl::Hit> hits;
        
         if(!igl::ray_mesh_intersect(s,dir,V,F,hits))
         {
             return false;
         }
         hit = hits[0];
         return true;
     };

     if(!shoot_ray(start,dir,hit))
     {
         return false;
     }

     Eigen::Vector3f bc;
     bc << 1.0-hit.u-hit.v, hit.u, hit.v;
     unsigned int fid = hit.id;
  
     long c;
     bc.maxCoeff(&c);
     unsigned int vid = F(fid,c);

     for(unsigned int qi = 0; qi < V.rows(); qi++) {
         if((V.row(qi) - V.row(vid)).norm() < radius) {
             verts.push_back(std::make_pair(gi,qi));
         }
     }

    //just pick the nearest vertex
    // Eigen::Vector3f mouse(win(0), win(1), win(2));
    // Eigen::Vector3f world_mouse;
    // igl::unproject(mouse,view,proj,viewport,world_mouse);
    // std::cout<<"world mouse:"<<world_mouse<<std::endl;
    
    // double dist=50.0;
    // int selected_idx = -1;
    // for(unsigned int qi = 0; qi < V.rows(); qi++) {
    //     double curr_dist = (V.row(qi).transpose() - world_mouse.cast<double>()).norm();
    //     if(curr_dist < dist) {
    //         selected_idx = qi;
    //         dist = curr_dist;
    //     }
    // }
    // std::cout<<"picked vertex:"<<V.row(selected_idx)<<std::endl;
    // if(selected_idx != -1){
    //     verts.push_back(selected_idx);
    // }

    return (verts.size() == 0 ? false : true);
}