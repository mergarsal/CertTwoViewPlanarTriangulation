#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

// triangulation solver
#include "TwoViewPlanarUtils.h"
// include types
#include "TwoViewPlanarTypes.h"
// include solver
#include "TwoViewPlanarDual.h"
// include SDP (Shor)
#include "TwoViewPlanarSDP.h"
// include certificate 
#include "TwoViewPlanarCert.h"


// point cloud 
#include "../utils/generatePointCloudPlanar.h"
#include "../utils/triangulationPlanar.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <chrono>  // timer


using namespace std;
using namespace Eigen;
using namespace TwoViewPlanar;
using namespace std::chrono;

size_t generateRandomIndex(size_t N)
{
        size_t n = N; 
        n = (((double)std::rand() - 0.5) / (double) RAND_MAX) * 2.0 * N;   
        return n; 
}

double distEuc(const Vector3 & X, const Vector3 & Y)
{
        double dist = 110.0; 
        dist = (X-Y).norm();
        return dist;
} 

double evalConstr(const Matrix3 & E, const Vector2 & p1, const Vector2 & p2)
{
        Vector3 o1, o2; 
        o1 << p1, 1; 
        o2 << p2, 1; 

        double res = o1.transpose() * E * o2; 
        return res;
}

int main(int argc, char** argv)
{
        std::cout << "Two view planar triangulation comparison\n"; 
        
        
        auto name_sol3d = "points_3D_example.txt"; 
        std::ofstream fsol3d(name_sol3d);                                         
                                         
                                         
        
          double noise = .5;         
          double max_parallax = 4.0;  // in meters
          double width_plane = 4;              // this is half the plane so: [-width, +width]
          double height_plane = 4;           
          size_t n_width = width_plane * 2;
          size_t n_height = height_plane * 2; 
          double focal_length = 512; 
          size_t size_img = 1024; 
          double d_plane = 1;

          size_t n_points = (n_width + 1) * (n_height + 1);
          
          
         std::srand(std::time(nullptr));
           
           
       Vector3 translation;
       Matrix3 rotation;
       
       
       Eigen::MatrixXd points_3D(3, n_points);
 
       // define struct with params
       PCPlanarParams str_in = PCPlanarParams();   
       str_in.focal_length = focal_length; 
       str_in.size_img = size_img;
       str_in.N_width = n_width; 
       str_in.N_height = n_height; 
       str_in.noise = noise; 
       str_in.width_plane = width_plane; 
       str_in.height_plane = height_plane;
       str_in.max_parallax = max_parallax; 
       str_in.d_plane = d_plane;
                                           
                                               
       // generate problem
       // You can set different parameters
       // str_in.max_angle = 0.0;                                                        
       // str_in.noise_trans = 0.0; 
       // str_in.noise_rot = 0.0; 
       // str_in.dir_parallax << 1, 0, 0; 
       PCPlanarRes str_out = generatePointCloudPlanar(str_in); // , generateTranslationStereo);     // generateTranslationSideways);  
       
                  
       // extract data
       translation = str_out.translation.normalized(); 
       rotation = str_out.rotation; 
       
       Matrix3 Tx; 
       Tx << 0, -translation(2), translation(1), translation(2), 0, -translation(0), -translation(1), translation(0), 0; 
       Matrix3 E = Tx * rotation;
       
       std::cout << "R:\n" << rotation << std::endl; 
       std::cout << "t:\n" << translation << std::endl; 
       
   
       
       Matrix3 H, Hp; 
                
                
       Vector3 n; 
       n << 0, 0, 1; 
       computeHomography(n, d_plane, str_out.R1, str_out.R2, str_out.t1, str_out.t2, focal_length, size_img, H, Hp);
                        

       
       
       Matrix3 K = Matrix3::Zero(); 
       K(0,0) = focal_length; 
       K(1,1) = focal_length; 
       K(0,2) = size_img / 2; 
       K(1,2) = size_img / 2; 
       K(2,2) = 1; 
       
       Matrix3 iK; 
       iK = K.inverse(); 
       Matrix3 F = iK.transpose() * E * iK;
       
       /* Save problem data */
       fsol3d << str_out.t1(0) << "," << str_out.t1(1) << "," << str_out.t1(2) << std::endl;              
       fsol3d << str_out.t2(0) << "," << str_out.t2(1) << "," << str_out.t2(2) << std::endl;   
           
       fsol3d << str_out.R1(0, 0) << "," << str_out.R1(0, 1)<< "," << str_out.R1(0, 2)<< std::endl;  
       fsol3d << str_out.R1(1, 0) << "," << str_out.R1(1, 1)<< "," << str_out.R1(1, 2)<< std::endl;  
       fsol3d << str_out.R1(2, 0) << "," << str_out.R1(2, 1)<< "," << str_out.R1(2, 2)<< std::endl;  
       
       fsol3d << str_out.R2(0, 0) << "," << str_out.R2(0, 1)<< "," << str_out.R2(0, 2)<< std::endl;  
       fsol3d << str_out.R2(1, 0) << "," << str_out.R2(1, 1)<< "," << str_out.R2(1, 2)<< std::endl;  
       fsol3d << str_out.R2(2, 0) << "," << str_out.R2(2, 1)<< "," << str_out.R2(2, 2)<< std::endl;  
       
       fsol3d << translation(0) << "," << translation(1) << "," << translation(2) << std::endl;   
       fsol3d << rotation(0, 0) << "," << rotation(0, 1)<< "," << rotation(0, 2)<< std::endl;  
       fsol3d << rotation(1, 0) << "," << rotation(1, 1)<< "," << rotation(1, 2)<< std::endl;  
       fsol3d << rotation(2, 0) << "," << rotation(2, 1)<< "," << rotation(2, 2)<< std::endl;
       
       fsol3d << Hp(0, 0) << "," << Hp(0, 1)<< "," << Hp(0, 2)<< std::endl;  
       fsol3d << Hp(1, 0) << "," << Hp(1, 1)<< "," << Hp(1, 2)<< std::endl;  
       fsol3d << Hp(2, 0) << "," << Hp(2, 1)<< "," << Hp(2, 2)<< std::endl;
       
       
       
       
       
       for (int i=0; i < 1; i++)
       {
                Vector3 p1, p2; 
                p1 = str_out.obs1.col(i); 
                p2 = str_out.obs2.col(i); 
                Vector3 point_i = str_out.points_3D.col(i);       
  
                Matrix3 Px = Matrix3::Zero(); 
                Px << 0, -p2(2), p2(1), p2(2), 0, -p2(0), -p2(1), p2(0), 0; 
                
                fsol3d << p1(0) << "," << p1(1) << ",1" << std::endl;   
                fsol3d << p2(0) << "," << p2(1) << ",1" << std::endl;   
                 
                // triangulate 3D point
                Vector3 Xi; 
                double depth_1, depth_2; 
                double error_svd = triangulatePoint(Xi, depth_1, depth_2, K * str_out.R1, 
                                                K * str_out.R2, K * str_out.t1, K * str_out.t2, p1, p2);
        
                double error_plane = errorPointPlane(Xi, n, d_plane);
                
                
                  
                   
                // Run solver 
                TwoViewPlanarClass solver_dual(Hp); 
                
                // solve!
                auto start_time_dual = high_resolution_clock::now();
                TwoViewPlanarResult res_dual = solver_dual.correctMatches(p1, p2); 
                auto time_init_dual = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_dual);   
                
                // Print result
                solver_dual.printResult(res_dual);    
                
                // Plot results 
                Vector3 q1 = res_dual.update_p1; 
                Vector3 q2 = res_dual.update_p2; 
                
                Matrix3 Qx = Matrix3::Zero(); 
                Qx << 0, -q2(2), q2(1), q2(2), 0, -q2(0), -q2(1), q2(0), 0; 
                
                fsol3d << q1(0) << "," << q1(1) << ",1" << std::endl;   
                fsol3d << q2(0) << "," << q2(1) << ",1" << std::endl; 
                
                
                std::cout << "Homography constraint dual: " << (Qx * Hp * q1).norm() << std::endl;
                
                 
                // triangulate 3D point
                Vector3 Xi_dual; 
                double depth_1_dual, depth_2_dual; 
                double error_svd_dual = triangulatePoint(Xi_dual, depth_1_dual, depth_2_dual, K * str_out.R1, 
                                                K * str_out.R2, K * str_out.t1, K * str_out.t2, q1, q2);
        
                double error_plane_dual = errorPointPlane(Xi_dual, n, d_plane);
                
                                            
                
                // Run SDP (Shor)
                
                TwoViewPlanarSDP solver_sdp = TwoViewPlanarSDP(Hp, p1, p2, false); 
                
                // solve SDP
                auto start_time_sdp = high_resolution_clock::now();
                TwoViewPlanarSDPResult res_sdp = solver_sdp.getResult(false);
                auto time_init_sdp = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_sdp);   
                
                // extract solution 
                solver_sdp.getSolutionFromResult(res_sdp, false); 
                
                Vector3 q1_sdp, q2_sdp; 
                q1_sdp = p1; 
                q2_sdp = p2;
                q1_sdp.block<2,1>(0,0) += res_sdp.delta_p.block<2,1>(0,0); 
                q2_sdp.block<2,1>(0,0) += res_sdp.delta_p.block<2,1>(2,0);
                Vector3 Xi_sdp; 
                double depth_1_sdp, depth_2_sdp; 
                double error_svd_sdp = triangulatePoint(Xi_sdp, depth_1_sdp, depth_2_sdp, K * str_out.R1, 
                                                K * str_out.R2, K * str_out.t1, K * str_out.t2, q1_sdp, q2_sdp);
        
                double error_plane_sdp = errorPointPlane(Xi_sdp, n, d_plane);
                
                Vector4 w_dual; 
                w_dual << res_dual.delta_p1(0), res_dual.delta_p1(1), res_dual.delta_p2(0), res_dual.delta_p2(1);
                
                std::cout << "[SDP] Dual cost = " << res_sdp.d_opt << std::endl; 
                std::cout << "[SDP] Primal cost = " << res_sdp.f_opt_sol << std::endl;  
                std::cout << "[SDP] Dual gap = " << res_sdp.d_opt - res_sdp.f_opt_sol << std::endl; 
                
                
                
                // Run polynomial solver 
                Vector2 point1, point2; 
                point1 << p1(0), p1(1); 
                point2 << p2(0), p2(1); 
                Vector2 q1_poly_nh, q2_poly_nh; 
                
                
                auto start_time_poly = high_resolution_clock::now();
                int n_roots_poly = correctMatchesPoly (Hp, point1, point2, & q1_poly_nh, & q2_poly_nh); 
                auto time_init_poly = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_poly);   
                
                
                
                Vector3 q1_poly, q2_poly; 
                q1_poly << q1_poly_nh, 1; 
                q2_poly << q2_poly_nh, 1; 
                
                fsol3d << (float) q1_poly_nh(0) << "," << q1_poly_nh(1) << ",1" << std::endl; 
                fsol3d << (float) q2_poly_nh(0) << "," << q2_poly_nh(1) << ",1" << std::endl; 
                 
                 
                Vector3 Xi_poly; 
                double depth_1_poly, depth_2_poly; 
                double error_svd_poly = triangulatePoint(Xi_poly, depth_1_poly, depth_2_poly, K * str_out.R1, 
                                                K * str_out.R2, K * str_out.t1, K * str_out.t2, q1_poly, q2_poly);
        
                double error_plane_poly = errorPointPlane(Xi_poly, n, d_plane);
                
                
                Matrix3 Hx = Matrix3::Zero(); 
                Hx << 0, -q2_poly(2), q2_poly(1), q2_poly(2), 0, -q2_poly(0), -q2_poly(1), q2_poly(0), 0; 
                               
                Vector3 diff_p1 = q1_poly - p1; 
                Vector3 diff_p2 = q2_poly - p2;   
                double cost_poly = diff_p1.dot(diff_p1) + diff_p2.dot(diff_p2);
                Vector4 w_poly; 
                w_poly << diff_p1(0), diff_p1(1), diff_p2(0), diff_p2(1); 
                
           
                
                
                /* Check certifier */
                Vector4 sol_w; 
                sol_w = w_dual; 
                double threshold_min = 1e-10;
                TwoViewPlanarCertRes res_cert = certifySolution(Hp, p1, p2, sol_w, threshold_min);
                
                /* Check sufficient condition */
                TwoViewPlanarSuffRes res_suff = certifySuffSolution(Hp, p1, p2, sol_w);
                
                
                std::cout << "Multipliers cert:\n" << res_cert.mult << std::endl; 
                std::cout << "Mult Dual:\n" << res_dual.lag_mult_1 << std::endl << res_dual.lag_mult_2 << std::endl; 
                
                
                // Run Kanatani's solver               
                
                Vector2 q1_iter_nh, q2_iter_nh; 
                
               
                auto start_time_iter = high_resolution_clock::now();
                int n_iter_kanatani = correctMatchesKanatani(Hp, point1, point2, 1, & q1_iter_nh, & q2_iter_nh); 
                auto time_init_iter = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_iter); 
                
                
                Vector3 q1_iter, q2_iter; 
                q1_iter << q1_iter_nh, 1; 
                q2_iter << q2_iter_nh, 1; 
                
                
                Vector3 Xi_iter; 
                double depth_1_iter, depth_2_iter; 
                double error_svd_iter = triangulatePoint(Xi_iter, depth_1_iter, depth_2_iter, K * str_out.R1, 
                                                K * str_out.R2, K * str_out.t1, K * str_out.t2, q1_iter, q2_iter);
        
                double error_plane_iter = errorPointPlane(Xi_iter, n, d_plane);
                
                
                
                                
                // check epipolar constraint                 
                
                std::cout << "[DUAL] Epipolar constraint: " << q2.transpose() * F * q1 << std::endl; 
                std::cout << "[DUAL] Error for svd: " << error_svd_dual << std::endl; 
                std::cout << "[DUAL] Euclidean distance: " << distEuc(point_i, Xi_dual) << std::endl;                 
                std::cout << "[DUAL] Error for plane: " << error_plane_dual << std::endl;   
                
                std::cout << "[SDP] Epipolar constraint: " << q2_sdp.transpose() * F * q1_sdp << std::endl; 
                std::cout << "[SDP] Error for svd: " << error_svd_sdp << std::endl; 
                std::cout << "[SDP] Euclidean distance: " << distEuc(point_i, Xi_sdp) << std::endl;                 
                std::cout << "[SDP] Error for plane: " << error_plane_sdp << std::endl;  
                
                
                std::cout << "[POLY] Epipolar constraint: " << q2_poly.transpose() * F * q1_poly << std::endl; 
                std::cout << "[POLY] Error for svd: " << error_svd_poly << std::endl; 
                std::cout << "[POLY] Euclidean distance: " << distEuc(point_i, Xi_poly) << std::endl;                 
                std::cout << "[POLY] Error for plane: " << error_plane_poly << std::endl;  
  
                
                
                std::cout << "GT point:\n" << point_i << std::endl; 
                std::cout << "Point SVD:\n" << Xi << std::endl; 
                std::cout << "Point Dual:\n" << Xi_dual << std::endl; 
                std::cout << "Point SDP:\n" << Xi_sdp << std::endl; 
                std::cout << "Point POLY:\n" << Xi_poly << std::endl; 
                std::cout << "Point ITER:\n" << Xi_iter << std::endl; 
                
                               
                
                std::cout << "[DUAL] Cost: " << res_dual.f_opt << std::endl; 
                std::cout << "[SDP]  Cost: " << res_sdp.f_opt_sol << std::endl;
                std::cout << "[POLY] Cost: " << cost_poly << std::endl;
                
   
                                
                std::cout << "Diff DUAL-POLY: " << (w_poly - w_dual).norm() / cost_poly << std::endl; 
                
                
                
                std::cout << "Time SDP:  " << (double) time_init_sdp.count() << std::endl;
                std::cout << "Time DUAL: " << (double) time_init_dual.count() << std::endl;
                std::cout << "Time POLY: " << (double) time_init_poly.count() << std::endl;
                std::cout << "Time ITER: " << (double) time_init_iter.count() << std::endl;
                

             
                
       }
       
       fsol3d.close();   
       
      
  return 0;

}  // end of main fcn
