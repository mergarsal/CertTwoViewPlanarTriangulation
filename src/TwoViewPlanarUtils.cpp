// include types
#include "TwoViewPlanarTypes.h"

// include header
#include "TwoViewPlanarUtils.h" 

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>  // for the determinant

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>


#include <chrono>  // timer

using namespace std::chrono; 


namespace TwoViewPlanar
{   
        HDecSparseResult decomposeHSparse(Matrix3 & H)
        {
                HDecSparseResult res = HDecSparseResult(); 
                
                auto start_time = high_resolution_clock::now();
                
                Matrix2 H1 = Matrix2::Zero(); 
                H1(0, 0) = -H(2, 0); 
                H1(0, 1) = -H(2, 1); 
                
                
                Eigen::JacobiSVD<Matrix2> svd(H1, Eigen::ComputeFullU | Eigen::ComputeFullV);
                
                res.s = svd.singularValues()(0); 
                res.U = svd.matrixV(); 
                res.UU = svd.matrixU().transpose(); 
                res.h31 = H(2, 0);
                res.h32 = H(2, 1);   
               
                res.change_A1 =  svd.matrixU()(0, 0) * svd.matrixU()(1, 1);  
                
                auto time_init = duration_cast<nanoseconds>(high_resolution_clock::now() - start_time); 
                
                res.time_dec = (double) time_init.count(); 
                return res;       
        
        }        
        
        
        
        double triangulatePoint(Vector3 & X, double & depth_1, double & depth_2, 
                                const Matrix3 & R1, const Matrix3 & R2, 
                                const Vector3 & t1, const Vector3 & t2, 
                                const Vector3 & p1, const Vector3 & p2)
                                
                                
        {
                
                Matrix6 P = Matrix6::Zero(); 
                
                P.block<3,3>(0, 0) = R1; 
                P.block<3,1>(0, 3) = t1; 
                
                P.block<3,3>(3, 0) = R2; 
                P.block<3,1>(3, 3) = t2; 
                
                P.block<3,1>(0, 4) = p1; 
                P.block<3,1>(3, 5) = p2;
                
                
                // compute linear system 
                Eigen::JacobiSVD<Matrix6> svd(P, Eigen::ComputeFullU | Eigen::ComputeFullV); 
                
                double error_svd = (svd.singularValues())(5); 
                Matrix6 V = svd.matrixV();
                
                Vector4 Xp = V.block<4,1>(0, 5); 
                X << Xp(0) / Xp(3), Xp(1) / Xp(3), Xp(2) / Xp(3); 
                depth_1 = V(4, 5); 
                depth_2 = V(5, 5);              
              
                
                return error_svd;
        }
        

        double errorPointPlane(const Vector3 & X, const Vector3 & n, const double d)
        {
                return (n.dot(X) - d);                
        
        }

}  // end of namespace
