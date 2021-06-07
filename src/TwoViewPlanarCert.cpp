#include "TwoViewPlanarCert.h"

#include <Eigen/Dense>  // for linear least-squares & eigenvalues

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>


namespace TwoViewPlanar
{



TwoViewPlanarCertRes certifySolution(const Matrix3& H, const Vector3 & p1, const Vector3 & p2, const Vector4 & sol_w, const double threshold_min)
{ 
 
                Matrix3 T1 = Matrix3::Zero(), T2 = Matrix3::Zero();
                T1 << 0, 0, 0, 0, 0, -1, 0, 1, 0; 
                T2 << 0, 0, 1, 0, 0, 0, -1, 0, 0; 
                Matrix23 Pro; 
                Pro << 1, 0, 0, 0, 1, 0; 
                
                // we normalize 
                
                Matrix3 H1, H2;       
                
                Matrix2 Hs1, Hs2;                 
                
                      
                H1 = T1.transpose() * H; 
                H2 = T2.transpose() * H;                                 
                                
                Hs1 = Pro * H1 * Pro.transpose(); 
                Hs2 = Pro * H2 * Pro.transpose(); 
                                
                // compute errors
                double c1  = p2.dot(H1 * p1); 
                double c2  = p2.dot(H2 * p1);                
                
                
                Vector2 f1 = 0.5 * Pro * H1 * p1;
                Vector2 f2 = 0.5 * Pro * H2 * p1; 
                
                Vector2 k1 = 0.5 * Pro * H1.transpose() * p2; 
                Vector2 k2 = 0.5 * Pro * H2.transpose() * p2; 
                
                Vector2 w1 = sol_w.block<2,1>(0,0); 
                Vector2 w2 = sol_w.block<2,1>(2,0);
                 
                  
                            
                // NOTE: This matrix has always full rank              
                Matrix42 As = Matrix42::Zero(); 
                As.block<2,1>(0,0) = 0.5 * Hs1.transpose() * w2 + k1;  
                As.block<2,1>(2,0) = 0.5 * Hs1 * w1 + f1;   
                As.block<2,1>(0,1) = 0.5 * Hs2.transpose() * w2 + k2;  
                As.block<2,1>(2,1) = 0.5 * Hs2 * w1 + f2;    
                
                
                Vector2 mult = As.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(sol_w);
                
               
                
                double s = 0.5 * std::sqrt(H(2, 0) * H(2, 0) + H(2, 1) * H(2, 1))*std::sqrt(mult(0)*mult(0) + mult(1)*mult(1)); 
                
                double min_eig = 1 - s;  
                
                double d_mult = - mult(0) * c1 - mult(1) * c2 - (mult(0) * k1 + mult(1) * k2).dot(w1) - (mult(0) * f1 + mult(1) * f2).dot(w2);
                
                // return result 
                TwoViewPlanarCertRes res = TwoViewPlanarCertRes(); 
                res.d_mult = d_mult; 
                res.min_eig = min_eig; 
                res.mult = mult; 
                // res.Hessian = Hess; 
                res.As = As; 
                res.is_opt = (min_eig > threshold_min) ? true : false; 
                
                return res;
                          
}

TwoViewPlanarSuffRes certifySuffSolution(const Matrix3& H, const Vector3 & p1, const Vector3 & p2, const Vector4 & sol_w)
{
        double s = std::sqrt(H(2, 0)* H(2,0) + H(2, 1) * H(2, 1)); 
        
         Matrix3 T1 = Matrix3::Zero(), T2 = Matrix3::Zero();
                T1 << 0, 0, 0, 0, 0, -1, 0, 1, 0; 
                T2 << 0, 0, 1, 0, 0, 0, -1, 0, 0; 
                Matrix23 Pro; 
                Pro << 1, 0, 0, 0, 1, 0; 
                
                // we normalize 
                
                Matrix3 H1, H2;    
                      
                H1 = T1.transpose() * H; 
                H2 = T2.transpose() * H;                                 
                                
                Vector2 f1 = 0.5 * Pro * H1 * p1;
                Vector2 f2 = 0.5 * Pro * H2 * p1; 
                
                Vector2 k1 = 0.5 * Pro * H1.transpose() * p2; 
                Vector2 k2 = 0.5 * Pro * H2.transpose() * p2; 
                
                Matrix42 B = Matrix42::Zero(); 
                B.block<2,1>(0,0) = k1; 
                B.block<2,1>(2,0) = f1;
                B.block<2,1>(0,1) = k2; 
                B.block<2,1>(2,1) = f2;  
        
                Eigen::JacobiSVD<Matrix42> svd(B, Eigen::ComputeFullU | Eigen::ComputeFullV);
              
                
                double rat_bs = (svd.singularValues())(0) / s; 
                
                
                
                bool is_opt = (std::sqrt(sol_w.dot(sol_w)) < rat_bs ) ? true: false;
        
                TwoViewPlanarSuffRes res = TwoViewPlanarSuffRes(); 
                res.s = s; 
                res.sb = (svd.singularValues())(0); 
                res.is_opt = is_opt; 

                return res; 
}



}  // end namesapce 
