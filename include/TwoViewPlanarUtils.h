#pragma once

#include <Eigen/Core>

// include types
#include "TwoViewPlanarTypes.h"

namespace TwoViewPlanar
{
        
         /** Struct for result: 
                sparse approach 
        **/
        struct HDecSparseResult
        {
                Matrix2 U = Matrix2::Zero(); 
                Matrix2 UU = Matrix2::Zero(); 
                double change_A1 = 1.0;                
                double s = 0; 
                double h31 = 0; 
                double h32 = 0; 
                double time_dec = 100;
        
                HDecSparseResult(){}; 
        };  // end of HDecSparseResult
   
                        
        HDecSparseResult decomposeHSparse(Matrix3 & H);         
        
        
        double triangulatePoint(Vector3 & X, double & depth_1, double & depth_2, 
                                const Matrix3 & R1, const Matrix3 & R2, 
                                const Vector3 & t1, const Vector3 & t2, 
                                const Vector3 & p1, const Vector3 & p2);
        
        
        double errorPointPlane(const Vector3 & X, const Vector3 & n, const double d);

}  // end of namespace
