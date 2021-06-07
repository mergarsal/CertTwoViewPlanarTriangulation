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
        std::cout << "Two view planar triangulation basic\n"; 
        
        
          double noise = 5;
          size_t n_width = 2;
          size_t n_height = 2;          
          double max_parallax = 2.;  // in meters
          double width_plane = 2;              // this is half the plane so: [-width, +width]
          double height_plane = 2; 
          double focal_length = 600; 
          size_t size_img = 1024; 
          double d_plane = 3;

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
       PCPlanarRes str_out = generatePointCloudPlanar(str_in); 
       // extract data
       translation = str_out.translation.normalized(); 
       rotation = str_out.rotation; 
       
       Matrix3 Tx; 
       Tx << 0, -translation(2), translation(1), translation(2), 0, -translation(0), -translation(1), translation(0), 0; 
       Matrix3 E = Tx * rotation;
       
       
 
       /* Compute homography matrices: 
       H: normalized = R + t * n^T / d
       Hp: in pixels: K * H * pinv(K)
       */
       
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
       

       for (int i=0; i < 1; i++)
       {
                Vector3 p1, p2; 
                p1 = str_out.obs1.col(i); 
                p2 = str_out.obs2.col(i); 
                Vector3 point_i = str_out.points_3D.col(i); 
                               
                
                Matrix3 Px = Matrix3::Zero(); 
                Px << 0, -p2(2), p2(1), p2(2), 0, -p2(0), -p2(1), p2(0), 0; 
                
                // triangulate 3D point
                Vector3 Xi; 
                double depth_1, depth_2; 
                double error_svd = triangulatePoint(Xi, depth_1, depth_2, K * str_out.R1, 
                                                K * str_out.R2, K * str_out.t1, K * str_out.t2, p1, p2);
        
                double error_plane = errorPointPlane(Xi, n, d_plane);
                
                
                // check epipolar constraint 
                std::cout << "Observation #" << i << std::endl; 
                std::cout << "Epipolar constraint: " << p2.transpose() * F * p1 << std::endl; 
                std::cout << "Homography constraint: " << Px * Hp * p1 << std::endl; 
                std::cout << "Error for svd: " << error_svd << std::endl; 
                std::cout << "Euclidean distance: " << distEuc(point_i, Xi) << std::endl;                 
                std::cout << "Error for plane: " << error_plane << std::endl;     
                   
                   
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
                
                
                // triangulate 3D point
                Vector3 Xi_dual; 
                double depth_1_dual, depth_2_dual; 
                double error_svd_dual = triangulatePoint(Xi_dual, depth_1_dual, depth_2_dual, K * str_out.R1, 
                                                K * str_out.R2, K * str_out.t1, K * str_out.t2, q1, q2);
        
                double error_plane_dual = errorPointPlane(Xi_dual, n, d_plane);
                
  
                /* Shor's relaxation */
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
              
                
                std::cout << "[SDP] Dual cost = " << res_sdp.d_opt << std::endl; 
                std::cout << "[SDP] Primal cost = " << res_sdp.f_opt_sol << std::endl;  
                std::cout << "[SDP] Dual gap = " << res_sdp.d_opt - res_sdp.f_opt_sol << std::endl; 
                
                

                /* Check certifier */
                Vector4 sol_w; 
                sol_w << res_dual.delta_p1, res_dual.delta_p2; 
                double threshold_min = 1e-10;
                TwoViewPlanarCertRes res_cert = certifySolution(Hp, p1, p2, sol_w, threshold_min);
                
                /* Check sufficient condition */
                TwoViewPlanarSuffRes res_suff = certifySuffSolution(Hp, p1, p2, sol_w); 
                
                
                
                
                
                
                
                
                std::cout << "[DUAL] Epipolar constraint: " << q2.transpose() * F * q1 << std::endl; 
                std::cout << "[DUAL] Error for svd: " << error_svd_dual << std::endl; 
                std::cout << "[DUAL] Euclidean distance: " << distEuc(point_i, Xi_dual) << std::endl;                 
                std::cout << "[DUAL] Error for plane: " << error_plane_dual << std::endl;   
                
                
                std::cout << "[SDP] Epipolar constraint: " << q2_sdp.transpose() * F * q1_sdp << std::endl; 
                std::cout << "[SDP] Error for svd: " << error_svd_sdp << std::endl; 
                std::cout << "[SDP] Euclidean distance: " << distEuc(point_i, Xi_sdp) << std::endl;                 
                std::cout << "[SDP] Error for plane: " << error_plane_sdp << std::endl;  
                
                std::cout << "GT point:\n" << point_i << std::endl; 
                std::cout << "Point SVD:\n" << Xi << std::endl; 
                std::cout << "Point Dual:\n" << Xi_dual << std::endl; 
                std::cout << "Point SDP:\n" << Xi_sdp << std::endl; 
                
                
                std::cout << "Multipliers SDP:\n" << res_sdp.dual_point << std::endl; 
                std::cout << "Multipliers Dual:\n" << res_dual.lag_mult_1 << std::endl << res_dual.lag_mult_2 << std::endl; 
                std::cout << "Multipliers Cert:\n" << res_cert.mult << std::endl; 
                
                
                std::cout << "Primal cost from SDP = " << res_sdp.f_opt << std::endl; 
                std::cout << "Dual cost from SDP   = " << res_sdp.d_opt << std::endl; 
                std::cout << "Dual cost from dual  = " << res_dual.rho << std::endl;  
                std::cout << "Dual cost from cert  = " << res_cert.d_mult << std::endl;   
                
                
                std::cout << "Time SDP: " << (double) time_init_sdp.count() << std::endl;
                std::cout << "Time DUAL: " << (double) time_init_dual.count() << std::endl;
                
                   
       }
       
        

 
  return 0;

}  // end of main fcn
