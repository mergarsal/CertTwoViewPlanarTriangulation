# pragma once

// include types
#include "TwoViewPlanarTypes.h"



namespace TwoViewPlanar
{
        
        struct TwoViewPlanarCertRes
        {
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
                double d_mult = 1000; 
                double min_eig = -1000; 
                Vector2 mult = Vector2::Zero(); 
                Matrix4 Hessian = Matrix4::Zero(); 
                Matrix42 As = Matrix42::Zero(); 
                bool is_opt = false; 
                
                /* Default constructor */
                TwoViewPlanarCertRes(){};
        
        };  // end of struct
        
        
        struct TwoViewPlanarSuffRes
        {
                double s = 1; 
                double sb = 0; 
                bool is_opt = false; 
                
                /* Default constructor */
                TwoViewPlanarSuffRes(){};
        };  // end of struct
        
        /** Certifier **/
        TwoViewPlanarCertRes certifySolution(const Matrix3& H, const Vector3 & p1, 
                                             const Vector3 & p2, const Vector4 & sol_w, 
                                             const double threshold_min);

        /** Sufficcient condition **/
        TwoViewPlanarSuffRes certifySuffSolution(const Matrix3& H, const Vector3 & p1, 
                                                 const Vector3 & p2, const Vector4 & sol_w,
                                                 const double threshold_min);

}  // end of namespace
