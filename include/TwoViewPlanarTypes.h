#pragma once


#include <Eigen/Core>

namespace TwoViewPlanar
{
                typedef Eigen::Matrix<double, 2, 1> Vector2;
                typedef Eigen::Matrix<double, 3, 1> Vector3;
                typedef Eigen::Matrix<double, 4, 1> Vector4;
                typedef Eigen::Matrix<double, 5, 1> Vector5;
                typedef Eigen::Matrix<double, 6, 1> Vector6;
               
                typedef Eigen::Matrix<double, 19, 1> Vector19;
                typedef Eigen::Matrix<double, 23, 1> Vector23;
                typedef Eigen::Matrix<double, 24, 1> Vector24;
                
                
                
                
                typedef Eigen::Matrix<double, 2, 2> Matrix2;
                typedef Eigen::Matrix<double, 3, 3> Matrix3;                
                typedef Eigen::Matrix<double, 4, 4> Matrix4;
                typedef Eigen::Matrix<double, 5, 5> Matrix5;
                typedef Eigen::Matrix<double, 6, 6> Matrix6;
                
                typedef Eigen::Matrix<double, 2, 3> Matrix23;
                typedef Eigen::Matrix<double, 4, 6> Matrix46;
                typedef Eigen::Matrix<double, 4, 2> Matrix42;
                



                /* lightweight struct for result */
                struct TwoViewPlanarResult
                {
                        double lag_mult_1 = 0; 
                        double lag_mult_2 = 0;
                        double rho = 0; 
                        
                        double f_opt = 0; 
                        
                        Vector3 update_p1; 
                        
                        Vector3 update_p2; 
                        
                        Vector2 delta_p1; 
                        
                        Vector2 delta_p2; 
                        
                        double time_init = 0.0; 
                        
                        double time_zero_coeff = 0.0;
                        
                        double time_zero = 0.0; 
                        
                        double time_recovery = 0.0; 
                        
                        double time_dec_H = 0.0; 
                        
                        int n_iter_zero = 0; 
                        
                        bool is_opt = false; 
                        
                        TwoViewPlanarResult(){};        
                
                };

}  // end of namespace
