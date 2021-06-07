// include types
#include "TwoViewPlanarTypes.h"

// include header
#include "TwoViewPlanarUtils.h" 

// include dual headers
#include "TwoViewPlanarDual.h"


#include <Eigen/Core>
#include <Eigen/SVD>


#include <vector>
#include <string>
#include <iostream>
#include <iomanip>


#include <chrono>  // timer
using namespace std::chrono;





namespace TwoViewPlanar
{   

        TwoViewPlanarResult TwoViewPlanarClass::correctMatches(const Vector3 & p1, const Vector3 & p2, bool use_affine)
        {
                
                // 0. Compute params 
                double c1 = 0;            // value of the error for the data p1, p2
                double c2 = 0; 
                Vector4 d1, d2; 
                double l1, l2;  // optimal values of the multipliers
                double c1_opt, c2_opt;    // value of the constraints for opt
                int n_iter_zero = 0;      // number iterations
                double a = s_ / 2; 
                
                Matrix3 T1 = Matrix3::Zero(), T2 = Matrix3::Zero();
                T1 << 0, 0, 0, 0, 0, -1, 0, 1, 0; 
                T2 << 0, 0, 1, 0, 0, 0, -1, 0, 0; 
                Matrix23 Pro; 
                Pro << 1, 0, 0, 0, 1, 0; 
                
                Matrix3 H1, H2;                
                
                auto start_time_init = high_resolution_clock::now();
                
                
                // assure the signs are like in the paper
                H1 = change_A1_ * T1.transpose() * H_; 
                H2 = T2.transpose() * H_;                                 
                                
                // compute errors
                c1  = p2.dot(H1 * p1); 
                c2  = p2.dot(H2 * p1);
                
             
                
                
                d1.block<2,1>(0, 0) = 0.5 * UU_.transpose() * Pro * H1 * p1;
                d2.block<2,1>(0, 0) = 0.5 * UU_.transpose() * Pro * H2 * p1; 
                
                d1.block<2,1>(2, 0) = 0.5 * U_.transpose() * Pro * H1.transpose() * p2; 
                d2.block<2,1>(2, 0) = 0.5 * U_.transpose() * Pro * H2.transpose() * p2; 
                
                   
                // Estimate initial guesses for multipliers
                double l_init_1 = 0;  
                double l_init_2 = 0; 
                
                if (use_affine)
                {
                        double dd1 = d1.dot(d1); 
                        double dd2 = d2.dot(d2); 
                        double dd12 = d1.dot(d2); 
                        
                        double det_dd = dd1*dd2 - dd12*dd12; 
                        l_init_1 = -0.5 * (dd2 * c1 - dd12 * c2) / det_dd; 
                        l_init_2 = -0.5 * (dd1 * c2 - dd12 * c1) / det_dd; 
                
                }
                else
                {
                        l_init_1 = -0.5 * c1 / (d1.dot(d1));  
                        l_init_2 = -0.5 * c2 / (d2.dot(d2)); 
                }
                
                
                
                auto time_init = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_init);   
                
                
                // 1. Compute parameters for Newton's
                
                auto start_time_init_newton = high_resolution_clock::now();
                // some useful variables
                double a2 = a  * a;
                double a3 = a2 * a; 
                double a4 = a3 * a; 
                double a5 = a4 * a; 
                double a6 = a5 * a; 
                
                
                // Simplification 
                // d1e0 = 0 = d2e1 = 0
                // d1e1 = - d2e0
                
                double d1e1 = d1(1), d1e2 = d1(2), d1e3 = d1(3); 
                double d2e2 = d2(2), d2e3 = d2(3); 
                
               
                
                
                // std::cout << "[DUAL] Computing coefficients for Newton's\n";
                /* Constraint 1* / 
                /**
                const1 = 
                ((2*a4*d1e0*d1e0 + 2*a4*d1e3*d1e3)*l1^5 + 
                (2*a4*d1e0*d1e1 + 2*a4*d1e0*d2e0 + 2*a4*d1e3*d2e3)*l1^4*l2 + 
                (c1*a4 + 2*d1e1*d1e2*a3)*l1^4 + 
                (4*a4*d1e0*d1e0 + 4*a4*d1e3*d1e3)*l1^3*l2^2 + 
                (- 4*a2*d1e0*d1e0 - 4*a2*d1e3*d1e3)*l1^3 +                 
                (6*a4*d1e0*d1e1 + 6*a4*d1e0*d2e0 - 2*a4*d1e1*d2e1 + 4*a4*d1e3*d2e3 - 2*a4*d2e0*d2e1)*l1^2*l2^3 +                 
                (2*a4*c1 + 6*a3*d1e1*d1e2 + 2*a3*d1e0*d2e2 + 2*a3*d1e2*d2e0 - 2*a3*d2e1*d2e2)*l1^2*l2^2 +                 
                (2*a2*d1e1*d2e1 - 4*a2*d1e0*d2e0 - 6*a2*d1e0*d1e1 + 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3)*l1^2*l2 + 
                (- 2*c1*a2 - 6*d1e1*d1e2*a)*l1^2 + 
                (2*a4*d1e1*d1e1 + 4*a4*d1e1*d2e0 + 2*a4*d1e3*d1e3 + 2*a4*d2e0*d2e0 - 2*a4*d2e1*d2e1 + 4*d1e0*a4*d2e1)*l1*l2^4 + 
                (4*a3*d1e1*d2e2 - 4*a3*d1e0*d1e2 + 4*a3*d1e2*d2e1 + 4*a3*d2e0*d2e2)*l1*l2^3 + 
                (- 2*a2*d1e0*d1e0 - 4*a2*d1e0*d2e1 - 4*a2*d1e1*d1e1 - 4*d2e0*a2*d1e1 - 2*a2*d1e2*d1e2 - 4*a2*d1e3*d1e3 + 2*a2*d2e1*d2e1 + 2*a2*d2e2*d2e2)*l1*l2^2 + 
                (4*a*d1e0*d1e2 - 4*a*d1e1*d2e2 - 4*a*d1e2*d2e1)*l1*l2 + (2*d1e0*d1e0 + 2*d1e1*d1e1 + 2*d1e2*d1e2 + 2*d1e3*d1e3)*l1 + 
                (2*a4*d1e1*d2e1 + 2*a4*d1e3*d2e3 + 2*a4*d2e0*d2e1)*l2^5 + 
                (a4*c1 - 2*a3*d1e0*d2e2 - 2*a3*d1e2*d2e0 + 2*a3*d2e1*d2e2)*l2^4 + 
                (- 2*a2*d1e0*d2e0 - 4*a2*d1e1*d2e1 - 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3 - 2*a2*d2e0*d2e1)*l2^3 + 
                (2*a*d1e0*d2e2 - 2*a2*c1 + 2*a*d1e2*d2e0 - 2*a*d2e1*d2e2)*l2^2 + 
                (2*d1e0*d2e0 + 2*d1e1*d2e1 + 2*d1e2*d2e2 + 2*d1e3*d2e3)*l2 + c1)
                /(a4*l1^4 + 2*a4*l1^2*l2^2 + a4*l2^4 - 2*a2*l1^2 - 2*a2*l2^2 + 1)

                **/
                // std::cout << "[DUAL] Going for C1\n";
                
                double c1l15 = 2*a4*d1e3*d1e3;                 
                double c1l14 = c1*a4 + 2*d1e1*d1e2*a3;                 
                double c1l13 = - 4*a2*d1e3*d1e3; 
                double c1l12 = - 2*c1*a2 - 6*d1e1*d1e2*a;
                double c1l1 = 2*d1e1*d1e1 + 2*d1e2*d1e2 + 2*d1e3*d1e3; 
                
                
                double c1l14l2 = 2*a4*d1e3*d2e3; 
                double c1l13l22 = 4*a4*d1e3*d1e3; 
                double c1l12l23 = 4*a4*d1e3*d2e3 ;
                double c1l12l22 = 2*a4*c1 + 4*a3*d1e1*d1e2;
               
                double c1l12l2 = 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3;                 
                double c1l1l24 = 2*a4*d1e1*d1e1 - 2*a4*d1e1*d1e1 + 2*a4*d1e3*d1e3;
                double c1l1l23 = 0; 
                double c1l1l22 = - 2*a2*d1e2*d1e2 - 4*a2*d1e3*d1e3 + 2*a2*d2e2*d2e2; 
                double c1l1l2 = - 4*a*d1e1*d2e2; 
                
                
                double c1l25 = 2*a4*d1e3*d2e3; 
                double c1l24 = a4*c1 + 2*a3*d1e2*d1e1; 
                double c1l23 = - 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3; 
                double c1l22 = - 2*a2*c1 - 2*a*d1e2*d1e1; 
                double c1l2 = 2*d1e2*d2e2 + 2*d1e3*d2e3; 
                
                // 
                Vector19 C1_coeff; 
                C1_coeff << c1l15, c1l14, c1l13, c1l12, c1l1, c1l25, c1l24, c1l23, c1l22, c1l2, c1l14l2, c1l13l22, c1l12l23, c1l12l22, c1l12l2, c1l1l24, c1l1l23, c1l1l22, c1l1l2; 
                
                /* Constraint 2 */ 
                /** 
                const2 = 
                ((2*a4*d1e0*d1e1 + 2*a4*d1e0*d2e0 + 2*a4*d1e3*d2e3)*l1^5 + 
                (- 2*a4*d1e0*d1e0 + 4*d2e1*a4*d1e0 + 2*a4*d1e1*d1e1 + 4*a4*d1e1*d2e0 + 2*a4*d2e0*d2e0 + 2*a4*d2e3*d2e3)*l1^4*l2 + 
                (a4*c2 - 2*a3*d1e0*d1e2 + 2*a3*d1e1*d2e2 + 2*a3*d1e2*d2e1)*l1^4 + 
                (6*a4*d1e1*d2e1 - 2*a4*d1e0*d2e0 - 2*a4*d1e0*d1e1 + 4*a4*d1e3*d2e3 + 6*a4*d2e0*d2e1)*l1^3*l2^2 + 
                (4*a3*d2e1*d2e2 - 4*a3*d1e0*d2e2 - 4*a3*d1e2*d2e0 - 4*a3*d1e1*d1e2)*l1^3*l2 + 
                (- 2*a2*d1e0*d1e1 - 4*a2*d1e0*d2e0 - 2*a2*d1e1*d2e1 - 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3)*l1^3 + 
                (4*a4*d2e1*d2e1 + 4*a4*d2e3*d2e3)*l1^2*l2^3 + 
                (2*a4*c2 + 2*a3*d1e0*d1e2 - 2*a3*d1e1*d2e2 - 2*a3*d1e2*d2e1 - 6*a3*d2e0*d2e2)*l1^2*l2^2 + 
                (2*a2*d1e0*d1e0 - 4*a2*d1e0*d2e1 + 2*a2*d1e2*d1e2 - 4*a2*d2e0*d2e0 - 4*d1e1*a2*d2e0 - 2*a2*d2e1*d2e1 - 2*a2*d2e2*d2e2 - 4*a2*d2e3*d2e3)*l1^2*l2 + 
                (2*a*d1e0*d1e2 - 2*a2*c2 - 2*a*d1e1*d2e2 - 2*a*d1e2*d2e1)*l1^2 + 
                (2*a4*d1e1*d2e1 + 2*a4*d1e3*d2e3 + 2*a4*d2e0*d2e1)*l1*l2^4 + 
                (2*a2*d1e0*d2e0 - 4*a2*d1e1*d2e1 + 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3 - 6*a2*d2e0*d2e1)*l1*l2^2 + 
                (4*a*d1e0*d2e2 + 4*a*d1e2*d2e0 - 4*a*d2e1*d2e2)*l1*l2 + 
                (2*d1e0*d2e0 + 2*d1e1*d2e1 + 2*d1e2*d2e2 + 2*d1e3*d2e3)*l1 + 
                (2*a4*d2e1*d2e1 + 2*a4*d2e3*d2e3)*l2^5 + 
                (a4*c2 - 2*a3*d2e0*d2e2)*l2^4 + 
                (- 4*a2*d2e1*d2e1 - 4*a2*d2e3*d2e3)*l2^3 + 
                (6*a*d2e0*d2e2 - 2*a2*c2)*l2^2 + 
                (2*d2e0*d2e0 + 2*d2e1*d2e1 + 2*d2e2*d2e2 + 2*d2e3*d2e3)*l2)
                **/
                                
                // std::cout << "[DUAL] Going for C2\n";
                
                double c2l15 = c1l14l2;                  
                double c2l14 = a4*c2 + 2*a3*d1e1*d2e2; 
                double c2l13 = - 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3;
                double c2l12 = - 2*a2*c2 - 2*a*d1e1*d2e2; 
                double c2l1 = c1l2;  
                
                double c2l14l2 = 2*a4*d1e1*d1e1 - 2*a4*d1e1*d1e1 + 2*a4*d2e3*d2e3;
                double c2l13l22 = 4*a4*d1e3*d2e3;                 
                double c2l13l2 = 0; 
                
                double c2l12l23 = 4*a4*d2e3*d2e3; 
                double c2l12l22 = 2*a4*c2 + 4*a3*d1e1*d2e2; 
                double c2l12l2 = 2*a2*d1e2*d1e2 - 2*a2*d2e2*d2e2 - 4*a2*d2e3*d2e3; 
                
                double c2l1l24 = c1l25; 
                double c2l1l22 = 2*a2*d1e2*d2e2 - 4*a2*d1e3*d2e3; 
                double c2l1l2 = -4*a*d1e2*d1e1; 
                
                
                double c2l25 = 2*a4*d2e3*d2e3; 
                double c2l24 = a4*c2 + 2*a3*d1e1*d2e2; 
                double c2l23 = - 4*a2*d2e3*d2e3; 
                double c2l22 = -6*a*d1e1*d2e2 - 2*a2*c2; 
                double c2l2 = 2*d1e1*d1e1 + 2*d2e2*d2e2 + 2*d2e3*d2e3; 
                
                Vector19 C2_coeff; 
                C2_coeff << c2l15, c2l14, c2l13, c2l12, c2l1, c2l25, c2l24, c2l23, c2l22, c2l2, c2l14l2, c2l13l22, c2l12l23, c2l12l22, c2l12l2, c2l1l24, c2l13l2, c2l1l22, c2l1l2; 
                
                
                /* Jacobian */
                /** J = [J11, J12; J21, J22]  **/ 
                
                /** J11 
                ((2*a6*d1e0*d1e0 + 2*a6*d1e3*d1e3)*l1^6 + 
                (6*a6*d1e0*d1e0 + 6*a6*d1e3*d1e3)*l1^4*l2^2 + 
                (- 6*a4*d1e0*d1e0 - 6*a4*d1e3*d1e3)*l1^4 + 
                (4*a6*d1e1*d2e1 - 4*a6*d1e0*d2e0 - 4*a6*d1e0*d1e1 + 4*a6*d2e0*d2e1)*l1^3*l2^3 + 
                (4*a5*d2e1*d2e2 - 4*a5*d1e0*d2e2 - 4*a5*d1e2*d2e0 - 4*a5*d1e1*d1e2)*l1^3*l2^2 + 
                (4*a4*d1e0*d1e1 - 4*a4*d1e1*d2e1 - 4*a4*d1e2*d2e2)*l1^3*l2 + 4*a3*d1e1*d1e2*l1^3 + 
                (12*a6*d1e0*d1e0 - 12*a6*d1e0*d2e1 - 6*a6*d1e1*d1e1 - 12*a6*d1e1*d2e0 + 6*a6*d1e3*d1e3 - 6*a6*d2e0*d2e0 + 6*a6*d2e1*d2e1)*l1^2*l2^4 + 
                (12*a5*d1e0*d1e2 - 12*a5*d1e1*d2e2 - 12*a5*d1e2*d2e1 - 12*a5*d2e0*d2e2)*l1^2*l2^3 + 
                (- 18*a4*d1e0*d1e0 + 12*a4*d1e0*d2e1 + 12*a4*d1e1*d1e1 + 12*d2e0*a4*d1e1 + 6*a4*d1e2*d1e2 - 12*a4*d1e3*d1e3 - 6*a4*d2e1*d2e1 - 6*a4*d2e2*d2e2)*l1^2*l2^2 + 
                (12*a3*d1e1*d2e2 - 12*a3*d1e0*d1e2 + 12*a3*d1e2*d2e1)*l1^2*l2 + 
                (6*a2*d1e0*d1e0 - 6*a2*d1e1*d1e1 - 6*a2*d1e2*d1e2 + 6*a2*d1e3*d1e3)*l1^2 + 
                (12*a6*d1e0*d1e1 + 12*a6*d1e0*d2e0 - 12*a6*d1e1*d2e1 - 12*a6*d2e0*d2e1)*l1*l2^5 + 
                (12*a5*d1e1*d1e2 + 12*a5*d1e0*d2e2 + 12*a5*d1e2*d2e0 - 12*a5*d2e1*d2e2)*l1*l2^4 + 
                (24*a4*d1e1*d2e1 - 12*a4*d1e0*d2e0 - 24*a4*d1e0*d1e1 + 12*a4*d1e2*d2e2 + 12*a4*d2e0*d2e1)*l1*l2^3 + 
                (12*a3*d2e1*d2e2 - 12*a3*d1e0*d2e2 - 12*a3*d1e2*d2e0 - 24*a3*d1e1*d1e2)*l1*l2^2 + 
                (12*a2*d1e0*d1e1 - 12*a2*d1e1*d2e1 - 12*a2*d1e2*d2e2)*l1*l2 + 
                12*a*d1e1*d1e2*l1 + 
                (2*a6*d1e1*d1e1 + 4*a6*d1e1*d2e0 + 2*a6*d1e3*d1e3 + 2*a6*d2e0*d2e0 - 2*a6*d2e1*d2e1 + 4*d1e0*a6*d2e1)*l2^6 + 
                (4*a5*d1e1*d2e2 - 4*a5*d1e0*d1e2 + 4*a5*d1e2*d2e1 + 4*a5*d2e0*d2e2)*l2^5 + 
                (- 2*a4*d1e0*d1e0 - 8*a4*d1e0*d2e1 - 6*a4*d1e1*d1e1 - 8*a4*d1e1*d2e0 - 2*a4*d1e2*d1e2 - 6*a4*d1e3*d1e3 - 2*a4*d2e0*d2e0 + 4*a4*d2e1*d2e1 + 2*a4*d2e2*d2e2)*l2^4 + 
                (8*a3*d1e0*d1e2 - 8*a3*d1e1*d2e2 - 8*a3*d1e2*d2e1 - 4*a3*d2e0*d2e2)*l2^3 + 
                (4*a2*d1e0*d1e0 + 4*a2*d1e0*d2e1 + 6*a2*d1e1*d1e1 + 4*d2e0*a2*d1e1 + 4*a2*d1e2*d1e2 + 6*a2*d1e3*d1e3 - 2*a2*d2e1*d2e1 - 2*a2*d2e2*d2e2)*l2^2 + 
                (4*a*d1e1*d2e2 - 4*a*d1e0*d1e2 + 4*a*d1e2*d2e1)*l2 - 2*d1e0*d1e0 - 2*d1e1*d1e1 - 2*d1e2*d1e2 - 2*d1e3*d1e3)
                **/
                
                // std::cout << "[DUAL] Going for J11\n";
                double j11l16 = 2*a6*d1e3*d1e3; 
                double j11l14 = - 6*a4*d1e3*d1e3; 
                double j11l13 = 4*a3*d1e1*d1e2; 
                double j11l12 = - 6*a2*d1e1*d1e1 - 6*a2*d1e2*d1e2 + 6*a2*d1e3*d1e3; 
                double j11l1 = 12*a*d1e1*d1e2; 
                
                 
                double j11l13l22 = 0; 
                double j11l13l2 = - 4*a4*d1e2*d2e2; 
                
                double j11l12l24 = 6*a6*d1e3*d1e3 ; 
                double j11l12l23 = 0; 
                double j11l12l22 = 6*a4*d1e2*d1e2 - 12*a4*d1e3*d1e3 - 6*a4*d2e2*d2e2; 
                double j11l12l2 = 12*a3*d1e1*d2e2 ; 
                double j11l14l22 =  6*a6*d1e3*d1e3; 
                double j11l13l23 = 0;
                double j11l1l25 = 0; 
                double j11l1l24 = 0; 
                double j11l1l23 = 12*a4*d1e2*d2e2; 
                double j11l1l22 =- 12*a3*d1e1*d1e2; 
                double j11l1l2 = - 12*a2*d1e2*d2e2; 
                
                
                double j11l26 = 2*a6*d1e3*d1e3; 
                double j11l25 = 0; 
                double j11l24 = - 2*a4*d1e2*d1e2 - 6*a4*d1e3*d1e3 + 2*a4*d2e2*d2e2; 
                double j11l23 = - 4*a3*d1e1*d2e2; 
                double j11l22 = 2 *a2*d1e1*d1e1  + 4*a2*d1e2*d1e2 + 6*a2*d1e3*d1e3 - 2*a2*d2e2*d2e2; 
                double j11l2 = 4*a*d1e1*d2e2; 
                
                
                Vector24 J11_coeff; 
                J11_coeff << j11l16,j11l14,j11l13,j11l12,j11l1, j11l26,j11l25,j11l24,j11l23,j11l22,j11l2, j11l13l22, j11l13l2, j11l12l24,j11l12l23,j11l12l22,j11l12l2, j11l14l22,j11l13l23,j11l1l25,j11l1l24,j11l1l23,j11l1l22,j11l1l2; 
                
                
                /** J22 
                
                ((- 2*a6*d1e0*d1e0 + 4*d2e1*a6*d1e0 + 2*a6*d1e1*d1e1 + 4*a6*d1e1*d2e0 + 2*a6*d2e0*d2e0 + 2*a6*d2e3*d2e3)*l1^6 + 
                (12*a6*d1e1*d2e1 - 12*a6*d1e0*d2e0 - 12*a6*d1e0*d1e1 + 12*a6*d2e0*d2e1)*l1^5*l2 + 
                (4*a5*d2e1*d2e2 - 4*a5*d1e0*d2e2 - 4*a5*d1e2*d2e0 - 4*a5*d1e1*d1e2)*l1^5 + 
                (6*a6*d1e0*d1e0 - 12*a6*d1e0*d2e1 - 6*a6*d1e1*d1e1 - 12*a6*d1e1*d2e0 - 6*a6*d2e0*d2e0 + 12*a6*d2e1*d2e1 + 6*a6*d2e3*d2e3)*l1^4*l2^2 + 
                (12*a5*d1e0*d1e2 - 12*a5*d1e1*d2e2 - 12*a5*d1e2*d2e1 - 12*a5*d2e0*d2e2)*l1^4*l2 + 
                (4*a4*d1e0*d1e0 - 8*a4*d1e0*d2e1 - 2*a4*d1e1*d1e1 - 8*a4*d1e1*d2e0 + 2*a4*d1e2*d1e2 - 6*a4*d2e0*d2e0 - 2*a4*d2e1*d2e1 - 2*a4*d2e2*d2e2 - 6*a4*d2e3*d2e3)*l1^4 + 
                (4*a6*d1e0*d1e1 + 4*a6*d1e0*d2e0 - 4*a6*d1e1*d2e1 - 4*a6*d2e0*d2e1)*l1^3*l2^3 + 
                (12*a5*d1e1*d1e2 + 12*a5*d1e0*d2e2 + 12*a5*d1e2*d2e0 - 12*a5*d2e1*d2e2)*l1^3*l2^2 + 
                (12*a4*d1e0*d1e1 + 24*a4*d1e0*d2e0 - 12*a4*d1e1*d2e1 + 12*a4*d1e2*d2e2 - 24*a4*d2e0*d2e1)*l1^3*l2 + 
                (4*a3*d1e1*d1e2 + 8*a3*d1e0*d2e2 + 8*a3*d1e2*d2e0 - 8*a3*d2e1*d2e2)*l1^3 + 
                (6*a6*d2e1*d2e1 + 6*a6*d2e3*d2e3)*l1^2*l2^4 + 
                (4*a5*d1e1*d2e2 - 4*a5*d1e0*d1e2 + 4*a5*d1e2*d2e1 + 4*a5*d2e0*d2e2)*l1^2*l2^3 +
                 (- 6*a4*d1e0*d1e0 + 12*a4*d1e0*d2e1 - 6*a4*d1e2*d1e2 + 12*a4*d2e0*d2e0 + 12*d1e1*a4*d2e0 - 18*a4*d2e1*d2e1 + 6*a4*d2e2*d2e2 - 12*a4*d2e3*d2e3)*l1^2*l2^2 + 
                 (12*a3*d1e1*d2e2 - 12*a3*d1e0*d1e2 + 12*a3*d1e2*d2e1 + 24*a3*d2e0*d2e2)*l1^2*l2 + 
                 (- 2*a2*d1e0*d1e0 + 4*a2*d1e0*d2e1 - 2*a2*d1e2*d1e2 + 6*a2*d2e0*d2e0 + 4*d1e1*a2*d2e0 + 4*a2*d2e1*d2e1 + 4*a2*d2e2*d2e2 + 6*a2*d2e3*d2e3)*l1^2 + 
                 (4*a4*d2e0*d2e1 - 4*a4*d1e2*d2e2 - 4*a4*d1e0*d2e0)*l1*l2^3 + 
                 (12*a3*d2e1*d2e2 - 12*a3*d1e2*d2e0 - 12*a3*d1e0*d2e2)*l1*l2^2 + 
                 (12*a2*d2e0*d2e1 - 12*a2*d1e2*d2e2 - 12*a2*d1e0*d2e0)*l1*l2 + 
                 (4*a*d2e1*d2e2 - 4*a*d1e2*d2e0 - 4*a*d1e0*d2e2)*l1 + 
                 (2*a6*d2e1*d2e1 + 2*a6*d2e3*d2e3)*l2^6 + 
                 (- 6*a4*d2e1*d2e1 - 6*a4*d2e3*d2e3)*l2^4 + 
                 (-4*a3*d2e0*d2e2)*l2^3 + 
                 (- 6*a2*d2e0*d2e0 + 6*a2*d2e1*d2e1 - 6*a2*d2e2*d2e2 + 6*a2*d2e3*d2e3)*l2^2 + 
                 (-12*a*d2e0*d2e2)*l2  
                **/
                // std::cout << "[DUAL] Going for J22\n";
                double j22l16 = 2*a6*d2e3*d2e3;                 
                double j22l15 = j11l13l22;   
                double j22l14 = 2*a4*d1e2*d1e2 - 2*a4*d2e2*d2e2 - 6*a4*d2e3*d2e3; 
                double j22l13 = -4*a3*d1e1*d1e2; 
                double j22l12 = 4*a2*d2e2*d2e2 + 6*a2*d2e3*d2e3;  
                double j22l1 = 4*a*d1e2*d1e1 ;    
                
                
                double j22l14l22 = 6*a6*d2e3*d2e3; 
                double j22l14l2 = j11l12l23; 
                double j22l15l2 = 0; 
                double j22l13l23 = 0; 
                double j22l13l22 = j11l1l24; 
                double j22l13l2 = 12*a4*d1e2*d2e2; 
                
                double j22l12l24 = 6*a6*d2e3*d2e3; 
                double j22l12l23 = j11l25; 
                double j22l12l22 = - 6*a4*d1e2*d1e2 + 6*a4*d2e2*d2e2 - 12*a4*d2e3*d2e3; 
                double j22l12l2 = - 12*a3*d1e1*d2e2; 
                double j22l1l23 = - 4*a4*d1e2*d2e2; 
                double j22l1l22 = 12*a3*d1e2*d1e1; 
                double j22l1l2 = - 12*a2*d1e2*d2e2; 
                 
                
                double j22l26 = 2*a6*d2e3*d2e3; 
                double j22l24 = - 6*a4*d2e3*d2e3; 
                double j22l23 = 4*a3*d1e1*d2e2; 
                double j22l22 = - 6*a2*d1e1*d1e1 - 6*a2*d2e2*d2e2 + 6*a2*d2e3*d2e3; 
                double j22l2 = 12*a*d1e1*d2e2; 
                
                
                Vector24 J22_coeff; 
                J22_coeff << j22l16, j22l15,j22l14,j22l13,j22l12,j22l1,j22l26,j22l24,j22l23,j22l22,j22l2, j22l14l22,j22l14l2,j22l15l2,j22l13l23,j22l13l22,j22l13l2,j22l12l24,j22l12l23,j22l12l22,j22l12l2,j22l1l23,j22l1l22,j22l1l2; 
                
                
                /** Entry J1,2
                
                ((2*a6*d1e0*d1e1 + 2*a6*d1e0*d2e0 + 2*a6*d1e3*d2e3)*l1^6 + 
                (12*a6*d1e0*d1e1 + 12*a6*d1e0*d2e0 - 6*a6*d1e1*d2e1 + 6*a6*d1e3*d2e3 - 6*a6*d2e0*d2e1)*l1^4*l2^2 + 
                (4*a5*d1e1*d1e2 + 4*a5*d1e0*d2e2 + 4*a5*d1e2*d2e0 - 4*a5*d2e1*d2e2)*l1^4*l2 + 
                (2*a4*d1e1*d2e1 - 6*a4*d1e0*d2e0 - 8*a4*d1e0*d1e1 + 2*a4*d1e2*d2e2 - 6*a4*d1e3*d2e3)*l1^4 + 
                (- 8*a6*d1e0*d1e0 + 16*a6*d1e0*d2e1 + 8*a6*d1e1*d1e1 + 16*a6*d1e1*d2e0 + 8*a6*d2e0*d2e0 - 8*a6*d2e1*d2e1)*l1^3*l2^3 + 
                (12*a5*d1e1*d2e2 - 12*a5*d1e0*d1e2 + 12*a5*d1e2*d2e1 + 12*a5*d2e0*d2e2)*l1^3*l2^2 + 
                (4*a4*d1e0*d1e0 - 8*a4*d1e0*d2e1 - 8*a4*d1e1*d1e1 - 8*d2e0*a4*d1e1 - 4*a4*d1e2*d1e2 + 4*a4*d2e1*d2e1 + 4*a4*d2e2*d2e2)*l1^3*l2 + 
                (4*a3*d1e0*d1e2 - 4*a3*d1e1*d2e2 - 4*a3*d1e2*d2e1)*l1^3 + 
                (12*a6*d1e1*d2e1 - 6*a6*d1e0*d2e0 - 6*a6*d1e0*d1e1 + 6*a6*d1e3*d2e3 + 12*a6*d2e0*d2e1)*l1^2*l2^4 + 
                (12*a5*d2e1*d2e2 - 12*a5*d1e0*d2e2 - 12*a5*d1e2*d2e0 - 12*a5*d1e1*d1e2)*l1^2*l2^3 + 
                (- 12*a4*d1e0*d2e0 - 12*a4*d1e1*d2e1 - 12*a4*d1e2*d2e2 - 12*a4*d1e3*d2e3)*l1^2*l2^2 + 12*a3*d1e1*d1e2*l1^2*l2 + 
                (6*a2*d1e0*d1e1 + 6*a2*d1e0*d2e0 + 6*a2*d1e3*d2e3)*l1^2 + 
                (4*a5*d1e0*d1e2 - 4*a5*d1e1*d2e2 - 4*a5*d1e2*d2e1 - 4*a5*d2e0*d2e2)*l1*l2^4 + 
                (4*a4*d1e0*d1e0 - 8*a4*d1e0*d2e1 + 4*a4*d1e2*d1e2 - 8*a4*d2e0*d2e0 - 8*d1e1*a4*d2e0 + 4*a4*d2e1*d2e1 - 4*a4*d2e2*d2e2)*l1*l2^3 + 
                (-12*a3*d2e0*d2e2)*l1*l2^2 + 
                (- 4*a2*d1e0*d1e0 + 8*a2*d1e0*d2e1 - 4*a2*d1e2*d1e2 - 4*a2*d2e1*d2e1 - 4*a2*d2e2*d2e2 + 8*d1e1*d2e0*a2)*l1*l2 + 
                (4*a*d1e1*d2e2 - 4*a*d1e0*d1e2 + 4*a*d1e2*d2e1)*l1 + 
                (2*a6*d1e1*d2e1 + 2*a6*d1e3*d2e3 + 2*a6*d2e0*d2e1)*l2^6 + 
                (2*a4*d1e0*d2e0 - 6*a4*d1e1*d2e1 + 2*a4*d1e2*d2e2 - 6*a4*d1e3*d2e3 - 8*a4*d2e0*d2e1)*l2^4 + 
                (4*a3*d1e0*d2e2 + 4*a3*d1e2*d2e0 - 4*a3*d2e1*d2e2)*l2^3 + 
                (6*a2*d1e1*d2e1 + 6*a2*d1e3*d2e3 + 6*a2*d2e0*d2e1)*l2^2 + 
                (4*a*d2e1*d2e2 - 4*a*d1e2*d2e0 - 4*a*d1e0*d2e2)*l2 
                **/ 
                // std::cout << "[DUAL] Going for J12\n";
                double j12l16 = 2*a6*d1e3*d2e3; 
                double j12l14 =  2*a4*d1e2*d2e2 - 6*a4*d1e3*d2e3;
                double j12l12 = 6*a2*d1e3*d2e3; 
                double j12l13 = - 4*a3*d1e1*d2e2;  
                double j12l1 = j11l2; 
                
                
                double j12l14l2 = 0 ; 
                double j12l14l22 =  6*a6*d1e3*d2e3; 
                double j12l13l23 = 0 ; 
                double j12l13l22 = 0; 
                double j12l13l2 = - 4*a4*d1e2*d1e2 + 4*a4*d2e2*d2e2;                 
                double j12l12l24 =  6*a6*d1e3*d2e3; 
                double j12l12l23 = 0; 
                double j12l12l22 = - 12*a4*d1e2*d2e2 - 12*a4*d1e3*d2e3; 
                double j12l12l2 = 12*a3*d1e1*d1e2;                 
                double j12l1l24 = 0; 
                double j12l1l23 = + 4*a4*d1e2*d1e2 - 4*a4*d2e2*d2e2; 
                double j12l1l22 = 12*a3*d1e1*d2e2; 
                double j12l1l2 = - 4*a2*d1e2*d1e2 - 4*a2*d2e2*d2e2 - 8*d1e1*d1e1*a2; 
                
                
                double j12l26 =  2*a6*d1e3*d2e3 ; 
                double j12l24 =  2*a4*d1e2*d2e2 - 6*a4*d1e3*d2e3; 
                double j12l23 = -4*a3*d1e2*d1e1 ; 
                double j12l22 = 6*a2*d1e3*d2e3 ; 
                double j12l2 = j22l1;           
                
                
                Vector23 J12_coeff; 
                J12_coeff << j12l16,j12l14,j12l13,j12l12,j12l1,j12l26,j12l24,j12l23,j12l22,j12l2, j12l14l2,j12l14l22,j12l13l23,j12l13l22,j12l13l2,j12l12l24,j12l12l23,j12l12l22,j12l12l2,j12l1l24,j12l1l23,j12l1l22,j12l1l2;
                
                
                
                
                auto time_init_newton = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_init_newton);   
                
               
                
                
               
                // Find zero in polynomial 
                auto start_time_zero = high_resolution_clock::now();                
                n_iter_zero = findZeroNewton(l_init_1, l_init_2, C1_coeff, C2_coeff, J11_coeff, J22_coeff, J12_coeff, c1, c2, a, d1.dot(d2), d1.dot(d1), d2.dot(d2), l1, l2, c1_opt, c2_opt);
                auto time_zero = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_zero);
                
                // 4. Recover solution from w
               
                
                auto start_time_rec = high_resolution_clock::now();
                
     
                Vector4 xp = Vector4::Zero();
                Vector4 wp = Vector4::Zero(); 
                // apply inverse
                double det_iH = a2*(l1*l1 + l2 * l2) - 1;
                xp(0, 0) = -(a2*l1*l1 -1) * (l2*d1e1) + (a2*l1*l2)*(l1*d1e1)     - (a*l2)*(l1*d1e2 + l2*d2e2);
                xp(1, 0) = (a2*l1*l2) * (-l2*d1e1)    + (a2*l2*l2 - 1)*(l1*d1e1) + (a*l1)*(l1*d1e2 + l2*d2e2);
                xp(2, 0) = (a*l2) * (l2*d1e1)         + (a*l1)*(l1*d1e1)        - (l1*d1e2 + l2*d2e2); 
                xp(0, 0) /= det_iH; 
                xp(1, 0) /= det_iH; 
                xp(2, 0) /= det_iH;
                xp(3, 0) = l1 * d1e3 + l2 * d2e3;
                     
                
                wp.block<2,1>(0, 0) = U_ * xp.block<2,1>(2, 0);
                wp.block<2,1>(2, 0) = UU_ * xp.block<2,1>(0, 0); 
                               
                // check optimality for dual               
                
                bool is_opt = (1 / (a2) >= l1*l1 + l2*l2) ? true : false;
               
                
                                 
                auto time_rec = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time_rec);
                
    
                // 4. Save result
                TwoViewPlanarResult res = TwoViewPlanarResult(); 
                res.lag_mult_1 = change_A1_ * l1;
                res.lag_mult_2 = l2; 
                res.rho = -l1*c1 - l2*c2 - (l1 * d1 + l2 * d2).dot(xp);
                res.f_opt = xp.dot(xp);  
                
                res.delta_p1 = wp.block<2, 1>(0, 0); 
                res.delta_p2 = wp.block<2, 1>(2, 0); 
                
                res.update_p1 = p1; 
                res.update_p1.block<2, 1>(0,0) += wp.block<2, 1>(0, 0); 
                res.update_p2 = p2; 
                res.update_p2.block<2, 1>(0,0) += wp.block<2, 1>(2, 0); 
        
                res.time_init = (double) time_init.count(); 
                res.time_zero = (double) time_zero.count(); 
                res.time_recovery = (double) time_rec.count(); 
                res.time_zero_coeff = (double) time_init_newton.count(); 
                res.time_dec_H = time_dec_H_;
                
                res.n_iter_zero = n_iter_zero; 
                
                res.is_opt = is_opt;  
                
                return res;           

        }
        
        
        void TwoViewPlanarClass::printResult(TwoViewPlanarResult & res)
        {
                std::cout << "[PRINT] Result from optimization\n"; 
                std::cout << "[PRINT] Lambda1 = " << res.lag_mult_1 << std::endl; 
                std::cout << "[PRINT] Lambda2 = " << res.lag_mult_2 << std::endl; 
                std::cout << "[PRINT] Rho = " << res.rho << std::endl; 
                std::cout << "[PRINT] Cost = " << res.f_opt << std::endl;
                std::cout << "[PRINT] Delta points.\nFor p1:\n" << res.delta_p1 << std::endl; 
                std::cout << "For p2:\n" << res.delta_p2 << std::endl; 
                
                std::cout << "[PRINT] Time init: " << res.time_init << std::endl; 
                std::cout << "[PRINT] Time zero init: " << res.time_zero_coeff << std::endl; 
                std::cout << "[PRINT] Time zero: " << res.time_zero << std::endl; 
                std::cout << "[PRINT] Time recovery: " << res.time_recovery << std::endl; 
                std::cout << "[PRINT] Time decomposition H: " << res.time_dec_H << std::endl; 
        
                std::cout << "[PRINT] Number iterations: " << res.n_iter_zero << std::endl; 
                std::cout << "[PRINT] Is optimal?: " << res.is_opt << std::endl; 
                std::cout << "[PRINT] Dual gap (f_opt - rho) = " << res.rho - res.f_opt << std::endl;
      
                return;
        }
        
        
        
        double TwoViewPlanarClass::findZeroNewton(double l_init_1, double l_init_2, 
                                                  Vector19 & C1_coeff, Vector19 &  C2_coeff, 
                                                  Vector24 & J11_coeff, Vector24 & J22_coeff, 
                                                  Vector23 & J12_coeff, 
                                                  double c1, double c2, 
                                                  double a, double dot_bh, 
                                                  double norm_b, double norm_h, 
                                                  double & l1, double & l2, 
                                                  double & c1_opt, double & c2_opt)
        {
        
                double l_up_1 = l_init_1, l_up_2 = l_init_2;
                double j11_up, j22_up, j12_up; 
                double i_j11, i_j22, i_j12; 
                
                
                // std::cout << "[ZERO] Original multiplier 1: " << l_up_1 << std::endl; 
                // std::cout << "[ZERO] Original multiplier 2: " << l_up_2 << std::endl; 
                        
                double f_up_1 = _evalConstr1(l_init_1, l_init_2, C1_coeff, c1, a);
                double f_up_2 = _evalConstr2(l_init_1, l_init_2, C2_coeff, c2, a); 
                // std::cout << "[ZERO] Constraint for init 1: " << f_up_1 << std::endl; 
                // std::cout << "[ZERO] Constraint for init 2: " << f_up_2 << std::endl;
         
                
                double inc_1 = 100; 
                double inc_2 = 100; 
                                
                bool stop_1 = false; 
                bool stop_2 = false; 
                
                double det_inv = 1; 
                
                int i = 0; 
                for (i=0; i < max_iter_; i++)                        
                {
                        // std::cout << "[ZERO] Iteration: " << i << std::endl; 
                        // 1. Compute gradient 
                        j11_up = _evalGrad11(l_up_1, l_up_2, J11_coeff, norm_b, a); 
                        j22_up = _evalGrad22(l_up_1, l_up_2, J22_coeff, norm_h, a); 
                        j12_up = _evalGrad12(l_up_1, l_up_2, J12_coeff, dot_bh, a); 
                        
                        // std::cout << "[ZERO] J11: " << j11_up << std::endl; 
                        // std::cout << "[ZERO] J22: " << j22_up << std::endl; 
                        // std::cout << "[ZERO] J12: " << j12_up << std::endl; 
                        // inverse 2x2 matrix
                        det_inv = j11_up * j22_up - j12_up * j12_up;  
                        
                        
                        i_j11 = j22_up / det_inv; 
                        i_j22 = j11_up / det_inv;
                        i_j12 = -j12_up / det_inv; 
                        
                        // 2. Update l
                        inc_1 = f_up_1 * i_j11 + f_up_2 * i_j12; 
                        inc_2 = f_up_1 * i_j12 + f_up_2 * i_j22;
                       
                       
                        
                        l_up_1 -= (inc_1);
                        l_up_2 -= (inc_2);  
                       
                        // std::cout << "[ZERO] New multiplier 1: " << l_up_1 << std::endl; 
                        // std::cout << "[ZERO] New multiplier 2: " << l_up_2 << std::endl; 
                         
                        // 3. Compute cost (also for next iteration) 
                        f_up_1 = _evalConstr1(l_up_1, l_up_2, C1_coeff, c1, a);
                        f_up_2 = _evalConstr2(l_up_1, l_up_2, C2_coeff, c2, a);  
                         
                        
                        // std::cout << "[ZERO] Constraint 1: " << f_up_1 << std::endl; 
                        // std::cout << "[ZERO] Constraint 2: " << f_up_2 << std::endl; 
                        
                        // std::cout << "[ZERO] Increment 1: " << inc_1 << std::endl; 
                        // std::cout << "[ZERO] Increment 2: " << inc_2 << std::endl; 
                        
                        
                        // 4. Check stop condition
                        stop_1 = ((f_up_1 < tol_) && (-tol_ < f_up_1)) || ( (inc_1 < tol_d_) && (-tol_d_ < inc_1) ); 
                        stop_2 = ((f_up_2 < tol_) && (-tol_ < f_up_2)) || ( (inc_2 < tol_d_) && (-tol_d_ < inc_2) ); 
                        
                        if ( stop_1 && stop_2)   break;                
                
                }
                
                // std::cout << "[ZERO] number of iterations: " << i + 1 << std::endl; 
                
                l1 = l_up_1;
                l2 = l_up_2; 
                c1_opt = f_up_1; 
                c2_opt = f_up_2;                
                
                return i; 
        }
        
        
        
        
        double TwoViewPlanarClass::_evalConstr1(double l1, double l2, Vector19 & coeff, double c1, double a)
        {
                /** 
                C1_coeff << c1l15, c1l14, c1l13, c1l12, c1l1, 
                [5] c1l25, c1l24, c1l23, c1l22, c1l2, 
                [10] c1l14l2, c1l13l22, c1l12l23, c1l12l22, c1l12l2, c1l1l24, c1l1l23, c1l1l22, c1l1l2; 
                **/
                double l12 = l1 *l1, l22 = l2  * l2; 
                double l13 = l12*l1, l23 = l22 * l2; 
                double l14 = l13*l1, l24 = l23 * l2; 
                double l15 = l14*l1, l25 = l24 * l2;
                double a2 = a* a; 
                double a4 = a2 * a2;
                                
                double den_c1 = a4*l14 + 2*a4*l12*l22 + a4*l24 - 2*a2*l12 - 2*a2*l22 + 1; 
                
                double c = coeff(0) * l15 + coeff(1)*l14 + coeff(2)*l13 + coeff(3)*l12 + coeff(4)*l1; 
                c += coeff(5) * l25 + coeff(6) * l24 + coeff(7) * l23 + coeff(8) * l22 + coeff(9) * l2; 
                c += coeff(10) * l14*l2 + coeff(11) * l13*l22 + coeff(12) * l12*l23 + coeff(13) * l12 * l22 + coeff(14) * l12*l2 + coeff(15)*l1*l24 + coeff(16)*l1*l23+ coeff(17)*l1*l22 + coeff(18)*l1*l2 + c1;                
        
        
                return (c / den_c1); 
        }



        double TwoViewPlanarClass::_evalConstr2(double l1, double l2, Vector19 & coeff, double c1, double a)
        {
                /** 
                C2_coeff << c2l15, c2l14, c2l13, c2l12, c2l1, 
                [5] c2l25, c2l24, c2l23, c2l22, c2l2, 
                [10] c2l14l2, c2l13l22, c2l12l23, c2l12l22, c2l12l2, c2l1l24, c2l13l2, c2l1l22, c2l1l2;  
                **/
                double l12 = l1 *l1, l22 = l2  * l2; 
                double l13 = l12*l1, l23 = l22 * l2; 
                double l14 = l13*l1, l24 = l23 * l2; 
                double l15 = l14*l1, l25 = l24 * l2;
                double a2 = a* a; 
                double a4 = a2 * a2;
                                
                double den_c1 = a4*l14 + 2*a4*l12*l22 + a4*l24 - 2*a2*l12 - 2*a2*l22 + 1; 
                
                double c = coeff(0) * l15 + coeff(1)*l14 + coeff(2)*l13 + coeff(3)*l12 + coeff(4)*l1; 
                c += coeff(5) * l25 + coeff(6) * l24 + coeff(7) * l23 + coeff(8) * l22 + coeff(9) * l2; 
                c += coeff(10) * l14*l2 + coeff(11) * l13*l22 + coeff(12) * l12*l23 + coeff(13) * l12 * l22 + coeff(14) * l12*l2 + coeff(15)*l1*l24 + coeff(16)*l13*l2+ coeff(17)*l1*l22 + coeff(18)*l1*l2 + c1;                
        
        
                return (c / den_c1); 
        }





        double TwoViewPlanarClass::_evalGrad11(double l1, double l2, Vector24 & coeff, double norm_d1_sq, double a)
                {
                        /** 
                        J11_coeff << j11l16,j11l14,j11l13,j11l12,j11l1, 
                        [5] j11l26,j11l25,j11l24,j11l23,j11l22,j11l2, 
                        [11] j11l13l22, j11l13l2, j11l12l24,j11l12l23,j11l12l22,j11l12l2, 
                        [17] j11l14l22,j11l13l23,j11l1l25,j11l1l24,j11l1l23,j11l1l22,j11l1l2;  
                        **/
                        
                        double l12 = l1 *l1, l22 = l2  * l2; 
                        double l13 = l12*l1, l23 = l22 * l2; 
                        double l14 = l13*l1, l24 = l23 * l2; 
                        double l15 = l14*l1, l25 = l24 * l2;
                        double l16 = l15*l1, l26 = l25 * l2;
                        
                        double a2 = a* a; 
                        double a4 = a2 * a2;
                        double a6 = a4 * a2;
                                        
                        double den_c1 = a6*l16 + 3*a6*l14*l22 + 3*a6*l12*l24 + a6*l26 - 3*a4*l14 - 6*a4*l12*l22 - 3*a4*l24 + 3*a2*l12 + 3*a2*l22 - 1; 
                        
                        double c = coeff(0) * l16 + coeff(1)*l14 + coeff(2)*l13 + coeff(3)*l12 + coeff(4)*l1; 
                        c += coeff(5) * l26 + coeff(6) * l25 + coeff(7) * l24 + coeff(8) * l23 + coeff(9) * l22 + coeff(10) * l2; 
                        c += coeff(11) * l13 * l22 + coeff(12) * l13*l2 + coeff(13) * l12*l24 + coeff(14) * l12*l23 + coeff(15) * l12*l22 + coeff(16) * l12*l2; 
                        c += coeff(17) * l14*l22 + coeff(18) * l13*l23 + coeff(19) * l1*l25 + coeff(20) * l1*l24 + coeff(21) * l1*l23 + coeff(22) * l1*l22 + coeff(23) * l1*l2 - 2 * norm_d1_sq;               
                                                
                        return (c / den_c1); 
                }
                
                
                
                
        
        double TwoViewPlanarClass::_evalGrad22(double l1, double l2, Vector24 & coeff, double norm_h1_sq, double a)
                {
                        /** 
                        J22_coeff << j22l16, j22l15,j22l14,j22l13,j22l12,j22l1,
                        [6] j22l26,j22l24,j22l23,j22l22,j22l2,     
                        [11] j22l14l22,j22l14l2,j22l15l2,j22l13l23,j22l13l22,j22l13l2,
                        [17] j22l12l24,j22l12l23,j22l12l22,j22l12l2,j22l1l23,j22l1l22,j22l1l2; 
                        **/
                        
                        double l12 = l1 *l1, l22 = l2  * l2; 
                        double l13 = l12*l1, l23 = l22 * l2; 
                        double l14 = l13*l1, l24 = l23 * l2; 
                        double l15 = l14*l1, l25 = l24 * l2;
                        double l16 = l15*l1, l26 = l25 * l2;
                        
                        double a2 = a* a; 
                        double a4 = a2 * a2;
                        double a6 = a4 * a2;
                                        
                        double den_c1 = a6*l16 + 3*a6*l14*l22 + 3*a6*l12*l24 + a6*l26 - 3*a4*l14 - 6*a4*l12*l22 - 3*a4*l24 + 3*a2*l12 + 3*a2*l22 - 1; 
                        
                        double c = coeff(0) * l16 + coeff(1)*l15 + coeff(2)*l14 + coeff(3)*l13 + coeff(4)*l12 + coeff(5) * l1; 
                        c += coeff(6) * l26 + coeff(7) * l24 + coeff(8) * l23 + coeff(9) * l22 + coeff(10) * l2; 
                        c += coeff(11) * l14 * l22 + coeff(12) * l14*l2 + coeff(13) * l15*l2 + coeff(14) * l13*l23 + coeff(15) * l13*l22 + coeff(16) * l13*l2; 
                        c += coeff(17) * l12*l24 + coeff(18) * l12*l23 + coeff(19) * l12*l22 + coeff(20) * l12*l2 + coeff(21) * l1*l23 + coeff(22) * l1*l22 + coeff(23) * l1*l2 - 2 * norm_h1_sq;               
                                                
                        return (c / den_c1); 
                }
                
                
                
                
         double TwoViewPlanarClass::_evalGrad12(double l1, double l2, Vector23 & coeff, double norm_bh, double a)
                {
                        /** 
                        J12_coeff << j12l16,j12l14,j12l13,j12l12,j12l1,
                        [5] j12l26,j12l4,j12l23,j12l22,j12l2,
                        [10] j12l14l2,j12l14l22,j12l13l23,j12l13l22,j12l13l2,j12l12l24,
                        [16] j12l12l23,j12l12l22,j12l12l2,j12l1l24,j12l1l23,j12l1l22,j12l1l2;
                        **/
                        
                        double l12 = l1 *l1, l22 = l2  * l2; 
                        double l13 = l12*l1, l23 = l22 * l2; 
                        double l14 = l13*l1, l24 = l23 * l2; 
                        double l15 = l14*l1, l25 = l24 * l2;
                        double l16 = l15*l1, l26 = l25 * l2;
                        
                        double a2 = a* a; 
                        double a4 = a2 * a2;
                        double a6 = a4 * a2;
                                        
                        double den_c1 = a6*l16 + 3*a6*l14*l22 + 3*a6*l12*l24 + a6*l26 - 3*a4*l14 - 6*a4*l12*l22 - 3*a4*l24 + 3*a2*l12 + 3*a2*l22 - 1; 
                        
                        
                        double c = coeff(0) * l16 + coeff(1)*l14 + coeff(2)*l13 + coeff(3)*l12 + coeff(4)*l1;  
                        c += coeff(5) * l26 + coeff(6) * l24 + coeff(7) * l23 + coeff(8) * l22 + coeff(9) * l2; 
                        c += coeff(10) * l14 * l2 + coeff(11) * l14*l22 + coeff(12) * l13*l23 + coeff(13) * l13*l22 + coeff(14) * l13*l2 + coeff(15) * l12*l24; 
                        c += coeff(16) * l12*l23 + coeff(17) * l12*l22 + coeff(18) * l12*l2 + coeff(19) * l1*l24 + coeff(20) * l1*l23 + coeff(21) * l1*l22 + coeff(22) * l1*l2 - 2 * norm_bh; 
                                      
                        return (c / den_c1); 
                }       
                

}   // end of namespace
