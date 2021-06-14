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
// include certificate 
#include "TwoViewPlanarCert.h"


// point cloud 
#include "../utils/generatePointCloudPlanar.h"
#include "../utils/triangulationPlanar.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <chrono>  // timer

// include helper 
#include "experimentsHelper.h"


using namespace std;
using namespace Eigen;
using namespace TwoViewPlanar;
using namespace std::chrono;



enum class methodGenPose{
        ORBITAL = 0, 
        LATERAL,
        STEREO, 
        FORWARD, 
        DIAGONAL
};  // end of enum class


Vector3 returnRandomTrans( double noise, const Vector3 & trans)
{       
        // This function returns the given translation
        return (trans);       
};


int main(int argc, char** argv)
{
          std::cout << "Generic test:\nTwo view planar triangulation comparison\n"; 
          
          
           /* Read params from input */ 
   
           string name_in_params = "basic_params.txt"; 
           
           /* Read the name of the file */
           if (argc > 1)
                name_in_params = argv[1]; 
                
          
          SceneOptions options; 
           
          std::cout << "Generic test file !\n";
          std::cout << "Input for the test: " << name_in_params << std::endl; 
          
          
          // read params from file
          bool valid_options = readOptionsFromFile(name_in_params, options);   
          
          // parameters for estimation
          /*
           double noise = 5.;
          size_t n_width = 2;
          size_t n_height = 2;          
          double max_parallax = 2.;  // in meters
          double width_plane = 2;              // this is half the plane so: [-width, +width]
          double height_plane = 2; 
          double focal_length = 600; 
          size_t size_img = 1024; 
          double d_plane = 3;

          size_t n_points = (n_width + 1) * (n_height + 1);
          */ 
                             
          
          std::srand(std::time(nullptr));
                    
          
    for (size_t noise_id=0; noise_id < options.n_noise; noise_id++)
    {
        double noise_i = options.noise[noise_id];

        for (size_t size_id = 0; size_id < options.n_size_img; size_id++)
        {
                double size_img_i = options.size_imgs[size_id];
        
        
                for (size_t par_id = 0; par_id < options.n_parallax; par_id++)
                {
                        double par_i = options.max_parallax[par_id];    
        
        
                        for (size_t focal_id = 0; focal_id < options.n_focal; focal_id++)
                        {
                                double focal_i = options.focal_length[focal_id];
                                
                               for (size_t d_c_id = 0; d_c_id < options.n_dist_centers; d_c_id++)
                               {
                                        double dist_center_i = options.dist_centers[d_c_id]; 
                                         
                                         // This file saves all our resutls
                                         auto name_f_3d = "err_3D_noise_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream f3d(name_f_3d);
                                         
                                         auto name_f_l2 = "err_l2_noise_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fl2(name_f_l2);
         
                                         
                                         
                                         auto name_f_l1 = "err_l1_noise_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fl1(name_f_l1);
                                         
                                         auto name_f_linfty = "err_linfty_noise_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream flinfty(name_f_linfty);
                                         
                                         auto name_f_times = "times_noise_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream ftime(name_f_times);
                                         
                                         auto name_f_dd = "times_diag_noise_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream ft_dual(name_f_dd);
                                         
                                         
                                         auto name_diff = "diff_noise_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fdiff(name_diff);
                                         
                                         
                                         
                                         
                                         auto name_sol = "sols_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fsol(name_sol);
                                                                              
                                         
                                         auto name_mu = "mu_min_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fmumin(name_mu);
                                         
                                         auto name_suff = "suff_cond_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fcond(name_suff);
                                         
                                         auto name_sol3d = "points_3D_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fsol3d(name_sol3d);
                                         
                                         auto name_mult = "multipliers_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fmult(name_mult);
                                         
                                         
                                         auto name_constr = "H_constr_" + std::to_string(noise_i) 
                                                                     + "_size_" + std::to_string((int)size_img_i) 
                                                                     + "_par_" + std::to_string(par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_center_" + std::to_string(dist_center_i) 
                                                                     + ".txt";
                                         std::ofstream fconstr(name_constr); 
                                         
                                       for (size_t n_iter = 0; n_iter < options.max_iter; n_iter++)
                                       {
                                                   
                                               Vector3 translation;
                                               Matrix3 rotation;
                                               
                                               
                                               size_t n_width = 18;
                                               size_t n_height = 14;          
                                               double width_plane = 9;              // this is half the plane so: [-width, +width]
                                               double height_plane = 7;  
                                               size_t n_points = (n_width + 1) * (n_height + 1);
                                               
                                               double d_plane = (double) dist_center_i; 
                                               
                                               Eigen::MatrixXd points_3D(3, n_points);
                                               Eigen::MatrixXd obs1(2, n_points);
                                               Eigen::MatrixXd obs2(2, n_points);
                                               
                                               
                                               // define struct with params
                                               PCPlanarParams str_in = PCPlanarParams();   
                                               str_in.focal_length = focal_i; 
                                               str_in.size_img = size_img_i;
                                               str_in.N_width = n_width; 
                                               str_in.N_height = n_height; 
                                               str_in.noise = noise_i; 
                                               str_in.width_plane = width_plane; 
                                               str_in.height_plane = height_plane;
                                               str_in.max_parallax = par_i; 
                                               str_in.d_plane = d_plane;                                               
                                               str_in.max_angle = options.max_rotation; 
          
          

                                               // Select pose generation 
                                               // param: options.method_trans in {1, 2, 3, 4}
                                               methodGenPose m_gen_pose = static_cast<methodGenPose>(options.method_trans);

                                               
                                               // generate problem
                                               PCPlanarRes str_out = PCPlanarRes(); 
                                               std::cout << "Selection method for pose generation\n"; 
                                               
                                               switch (m_gen_pose)
                                               {
                                                case methodGenPose::ORBITAL:
                                                {
                                                        std::cout << "[GENERAL CAMERA]\n";
                                                        
                                                        str_out = generatePointCloudPlanar(str_in); // , returnRandomTrans, generateOrbitalRotation); 
                                                }
                                                break;
                                                
                                                case methodGenPose::LATERAL: 
                                                        std::cout << "[LATERAL CAMERA]\n";
                                                        str_in.max_angle = 0.0;                                                        
                                                        str_in.noise_trans = 0.01; 
                                                        str_in.noise_rot = 0.01; 
                                                        str_out = generatePointCloudPlanar(str_in, generateTranslationSideways);                                                 
                                                break;
                                                
                                                case methodGenPose::STEREO: 
                                                        std::cout << "[STEREO CAMERA]\n";
                                                        str_in.max_angle = 0.0;                                                        
                                                        str_in.noise_trans = 0.01; 
                                                        str_in.noise_rot = 0.01; 
                                                        str_out = generatePointCloudPlanar(str_in, generateTranslationStereo);                                                 
                                                break;
                                                
                                                case methodGenPose::FORWARD: 
                                                        std::cout << "[FORWARD CAMERA]\n";
                                                        str_in.max_angle = 0.0;                                                         
                                                        str_in.noise_trans = 0.01; 
                                                        str_in.noise_rot = 0.01;
                                                        str_out = generatePointCloudPlanar(str_in, generateTranslationForward);                                                 
                                                break;
        
                                                case methodGenPose::DIAGONAL: 
                                                        std::cout << "[DIAGONAL CAMERA]\n";
                                                        str_in.max_angle = 0.0; 
                                                        
                                                        str_in.noise_trans = 0.01; 
                                                        str_in.noise_rot = 0.01;
                                                        str_out = generatePointCloudPlanar(str_in, generateTranslationOblique);                                                 
                                                break;
        
                                                default:
                                                        std::cout << "[LATERAL CAMERA]\n";
                                                        str_in.max_angle = 0.0; 
                                                        str_in.noise_trans = 0.01; 
                                                        str_in.noise_rot = 0.01;
                                                        str_out = generatePointCloudPlanar(str_in, generateTranslationSideways);                                       
                                                break;
                                               
                                               }
                                                
                                               // extract data
                                               translation = str_out.translation.normalized(); 
                                               rotation = str_out.rotation; 
                                               
                                               std::cout << "Translation:\n" << translation << std::endl; 
                                               std::cout << "Rotation:\n" << rotation << std::endl; 
                                               
                                               Matrix3 Tx; 
                                               Tx << 0, -translation(2), translation(1), translation(2), 0, -translation(0), -translation(1), translation(0), 0; 
                                               Matrix3 E = Tx * rotation;
                                               
                                               
                                               Matrix3 H, Hp; 
                
                                               Vector3 n; 
                                               n << 0, 0, 1; 
                                               computeHomography(n, d_plane, str_out.R1, str_out.R2, str_out.t1, str_out.t2, focal_i, size_img_i, H, Hp);
                                               Matrix3 K = Matrix3::Zero(); 
                                               K(0,0) = focal_i; 
                                               K(1,1) = focal_i; 
                                               K(0,2) = size_img_i / 2; 
                                               K(1,2) = size_img_i / 2; 
                                               K(2,2) = 1; 
                                               
                                               Matrix3 iK; 
                                               iK = K.inverse(); 
                                               Matrix3 F = iK.transpose() * E * iK;                                                    
                                               
                                                                                              
                                               
                                               int idx = 0; 
                                               for (idx = 0; idx < n_points; idx++)
                                               {
                                                        
                                                        
                                                        Vector3 p1, p2;    // observations in sphere
                                                        Vector2 z1, z2;  // observations in Z plane (normalized)
                                                        Vector3 P;         // 3D point
                                                        
                                                        // 0. clean data
                                                        P.setZero(); 
                                                        
                                                        p1.setZero(); 
                                                        p2.setZero(); 
                                                        
                                                        z1.setZero(); 
                                                        z2.setZero(); 
                                                        // a. Extract observations
                                                        
                                                        p1 = str_out.obs1.col(idx);  
                                                        p2 = str_out.obs2.col(idx);   
                                                                     
                                                        P = str_out.points_3D.col(idx);
                                                        
                                                       
                                                        
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
                                                        
                                                     
                                                        
                                                        // Plot results 
                                                        Vector3 q1 = res_dual.update_p1; 
                                                        Vector3 q2 = res_dual.update_p2; 
                                                                                                                
                                                        Vector4 w_dual; 
                                                        w_dual << res_dual.delta_p1, res_dual.delta_p2; 
                                                        // triangulate 3D point
                                                        Vector3 Xi_dual; 
                                                        double depth_1_dual, depth_2_dual; 
                                                        double error_svd_dual = triangulatePoint(Xi_dual, depth_1_dual, depth_2_dual, K * str_out.R1, 
                                                                                        K * str_out.R2, K * str_out.t1, K * str_out.t2, q1, q2);
                                                
                                                        double error_plane_dual = errorPointPlane(Xi_dual, n, d_plane);
                
               
                                                        
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
                                                        
                                                        Vector4 w_poly; 
                                                        w_poly << q1_poly_nh - point1, q2_poly_nh - point2;
                                                         
                                                        Vector3 Xi_poly; 
                                                        double depth_1_poly, depth_2_poly; 
                                                        double error_svd_poly = triangulatePoint(Xi_poly, depth_1_poly, depth_2_poly, K * str_out.R1, 
                                                                                        K * str_out.R2, K * str_out.t1, K * str_out.t2, q1_poly, q2_poly);
                                                
                                                        double error_plane_poly = errorPointPlane(Xi_poly, n, d_plane);
                                                        
                
                                                          
                                                        /* Check certifier */
                                                        Vector4 sol_w; 
                                                        sol_w = w_dual; 
                                                        double threshold_min = 1e-10;
                                                        
                                                        auto start_time_cert = high_resolution_clock::now();
                                                        TwoViewPlanarCertRes res_cert = certifySolution(Hp, p1, p2, sol_w, threshold_min);
                                                        auto time_init_cert = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_cert);  
                                                        
                                                        
                                                        /* Check sufficient condition */
                                                        auto start_time_suff = high_resolution_clock::now();                                                        
                                                        TwoViewPlanarSuffRes res_suff = certifySuffSolution(Hp, p1, p2, sol_w, threshold_min);
                                                        auto time_init_suff = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_suff);  
                
                
                                                        
                                                        // Run Kanatani's solver
                                                        
                                                        
                                                        Vector2 q1_iter_nh, q2_iter_nh; 
                                                        
                                                        auto start_time_iter = high_resolution_clock::now();
                                                        int n_iter_kanatani = correctMatchesKanatani(Hp, point1, point2, 1, & q1_iter_nh, & q2_iter_nh); 
                                                        auto time_init_iter = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_iter); 
                                                        
                                                        
                                                        Vector3 q1_iter, q2_iter; 
                                                        q1_iter << q1_iter_nh, 1; 
                                                        q2_iter << q2_iter_nh, 1; 
                                                        
                                                        Vector4 w_iter; 
                                                        w_iter << q1_iter_nh - point1, q2_iter_nh - point2;
                                                        
                                                        Vector3 Xi_iter; 
                                                        double depth_1_iter, depth_2_iter; 
                                                        double error_svd_iter = triangulatePoint(Xi_iter, depth_1_iter, depth_2_iter, K * str_out.R1, 
                                                                                        K * str_out.R2, K * str_out.t1, K * str_out.t2, q1_iter, q2_iter);
                                                
                                                        double error_plane_iter = errorPointPlane(Xi_iter, n, d_plane);
                                                                                                        
                                                                                                              
                                        
                                                        /* Genera solver HS */
                                                        // project to normalized plane
                                                        
                                                        Vector2 hs_p1_2, hs_p2_2; 
                                                        auto start_t_hs = high_resolution_clock::now();
                                                        int n_roots = FindObsHS(F.transpose(), point1, point2, &hs_p1_2, &hs_p2_2);
                                                        auto time_hs = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_hs); 
                                                        
                                                        Vector4 w_hs; 
                                                        // NOTE: the difference is computed the other way in this solver
                                                        w_hs <<  (point1 - hs_p1_2), (point2 - hs_p2_2); 
                                                        Vector3 hs3_p1, hs3_p2; 
                                                        hs3_p1 << hs_p1_2, 1; 
                                                        hs3_p2 << hs_p2_2, 1; 
                                                        double depth_hs1, depth_hs2; 
                                                        Vector3 Xi_hs; 
                                                        double error_svd_hs = triangulatePoint(Xi_hs, depth_hs1, depth_hs2, K * str_out.R1, 
                                                                                        K * str_out.R2, K * str_out.t1, K * str_out.t2, hs3_p1, hs3_p2);
                                                                                        
                                                        
                
                
                                                                                                
                                                        
                                                        /* Save results */ 
                                                        // Euclidean distance 3D point
                                                        f3d << (Xi-P).norm() / par_i << ","; 
                                                        f3d << (Xi_dual-P).norm() / par_i << ","; 
                                                        f3d << (Xi_poly-P).norm() / par_i << ","; 
                                                        f3d << (Xi_iter-P).norm() / par_i << ",";                     
                                                        f3d << (Xi_hs-P).norm() / par_i; 
                                                        f3d << std::endl;
                                                        
                                                        
                                                        // l2 norm for observations 
                                                        fl2 << w_dual.norm()   << ","; 
                                                        fl2 << w_poly.norm()   << ",";  
                                                        fl2 << w_iter.norm()   << ",";                     
                                                        fl2 << w_hs.norm(); 
                                                        fl2 << std::endl;
                                                        
                                                        // L1 norm for observations 
                                                        fl1 << w_dual.lpNorm<1>()    << ","; 
                                                        fl1 << w_poly.lpNorm<1>()   << ",";  
                                                        fl1 << w_iter.lpNorm<1>()    << ",";                     
                                                        fl1 << w_hs.lpNorm<1>(); 
                                                        fl1 << std::endl;
                                                        
                                                        // Linfty norm for observations                                                         
                                                        flinfty << w_dual.lpNorm<Infinity>() << ","; 
                                                        flinfty << w_poly.lpNorm<Infinity>() << ","; 
                                                        flinfty << w_iter.lpNorm<Infinity>() << ",";                      
                                                        flinfty << w_hs.lpNorm<Infinity>(); 
                                                        flinfty << std::endl;
                                                        
                                                        /* Save diff wrt HS */
                                                        fdiff << (w_dual - w_poly).norm(); 
                                                        fdiff << std::endl; 
                                                       
                                                        
                                                        /** Minimum eigenvalue **/ 
                                                        fmumin << res_cert.min_eig << std::endl; 
                                                        
                                                        /* Suff. cond **/ 
                                                        
                                                        fcond << res_suff.sb << ","; 
                                                        fcond << res_suff.s; 
                                                        fcond << std::endl; 
                                                        
                                                        /** Times for all the solvers **/
                                                        ftime << (double) time_init_dual.count() << ",";                                                      
                                                        ftime << (double) time_init_poly.count() << ","; 
                                                        ftime << (double) time_init_iter.count() << ",";                                                         
                                                        ftime << (double) time_hs.count() << ","; 
                                                        ftime << (double) time_init_cert.count() << ","; 
                                                        ftime << (double) time_init_suff.count(); 
                                                        ftime << std::endl;        
                                                        
                                                        
                                                        
                                                        /** Times for dual solver **/
                                                        ft_dual <<  res_dual.time_dec_H << ","; 
                                                        ft_dual <<  res_dual.time_init << ","; 
                                                        ft_dual <<  res_dual.time_zero_coeff << ","; 
                                                        ft_dual <<  res_dual.time_zero << ","; 
                                                        ft_dual <<  res_dual.time_recovery << ","; 
                                                        ft_dual <<  res_dual.n_iter_zero << std::endl;        
                                                            
                                                         
                                                       
                                                      
                                                        
                                                        /** Solution **/
                                                        fsol << w_dual(0) << "," << w_dual(1) << "," << w_dual(2) << "," << w_dual(3) << std::endl; 
                                                        fsol << w_poly(0) << "," << w_poly(1) << "," << w_poly(2) << "," << w_poly(3) << std::endl; 
                                                        fsol << w_iter(0) << "," << w_iter(1) << "," << w_iter(2) << "," << w_iter(3) << std::endl; 
                                                        
                                                        
                                                        /**  Reconstructed point  **/
                                                        fsol3d << P(0) << "," << P(1) << "," << P(2) << std::endl; 
                                                        fsol3d << Xi(0) << "," << Xi(1) << "," << Xi(2) << std::endl; 
                                                        fsol3d << Xi_dual(0) << "," << Xi_dual(1) << "," << Xi_dual(2) << std::endl; 
                                                        fsol3d << Xi_poly(0) << "," << Xi_poly(1) << "," << Xi_poly(2) << std::endl; 
                                                        fsol3d << Xi_iter(0) << "," << Xi_iter(1) << "," << Xi_iter(2) << std::endl;                   
                                                        fsol3d << Xi_hs(0) << "," << Xi_hs(1) << "," << Xi_hs(2) << std::endl; 
                                                        
                                                        Vector2 mult_dual; 
                                                        mult_dual << res_dual.lag_mult_1, res_dual.lag_mult_2;
                                                        fmult << res_dual.lag_mult_1 << "," << res_dual.lag_mult_2; 
                                                        fmult << "," << res_cert.mult(0) << "," << res_cert.mult(1) << ","; 
                                                        fmult << (mult_dual - res_cert.mult).norm();
                                                        fmult << std::endl; 
                                                        
                                                                                                                /* Constraint */
                                                        Matrix3 Px = Matrix3::Zero(); 
                                                        /* data */
                                                        Px << 0, -p2(2), p2(1), p2(2), 0, -p2(0), -p2(1), p2(0), 0;                                                         
                                                        fconstr << (Px * Hp * p1).norm() << ","; 
                                                        /* dual */
                                                        Px << 0, -q2(2), q2(1), q2(2), 0, -q2(0), -q2(1), q2(0), 0; 
                                                        fconstr << (Px * Hp * q1).norm() << ","; 
                                                        /* poly */
                                                        Px << 0, -q2_poly(2), q2_poly(1), q2_poly(2), 0, -q2_poly(0), -q2_poly(1), q2_poly(0), 0; 
                                                        fconstr << (Px * Hp * q1_poly).norm() << ","; 
                                                        /* iter */
                                                        Px << 0, -q2_iter(2), q2_iter(1), q2_iter(2), 0, -q2_iter(0), -q2_iter(1), q2_iter(0), 0; 
                                                        fconstr << (Px * Hp * q1_iter).norm() << ","; 
                                                        /* poly */
                                                        Px << 0, -hs3_p2(2), hs3_p2(1), hs3_p2(2), 0, -hs3_p2(0), -hs3_p2(1), hs3_p2(0), 0; 
                                                        fconstr << (Px * Hp * hs3_p1).norm() << ",";  
                                                        fconstr << hs3_p2.transpose() * F * hs3_p1; 
                                                        fconstr << std::endl; 
                                                        
                                                        
                                                }  // end of each point             
                                                } // end of each iteration 
                                                
                                                f3d.close();      // error 3d point
                                                fl2.close();      // error l2 obs
                                                fl1.close();      // error l1 obs
                                                flinfty.close();  // error linfty obs                                                
                                                ftime.close();    // time 
                                                ft_dual.close();  // time dual diagonal                                      
                                                fdiff.close();    // diff wrt HS
                                                fsol.close();     // solutions L2
                                                fmumin.close();   // mu min
                                                fcond.close();    // suff. cond
                                                fsol3d.close();   // reconstructed 3D points
                                                fmult.close();    // multipliers
                                                fconstr.close();   // constraint
                                                
                                }  // end for dist center
      
                         }  // enf for focal
       
                 }  // end for parallax
      
        }  // end for image size
      
      }  // end for noise              
       
  return 0;

}  // end of main fcn
