#pragma once

#include <eigen3/Eigen/Dense>

#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

// include types
#include "TwoViewPlanarTypes.h"
// include basic fcns
#include "TwoViewPlanarUtils.h"



namespace TwoViewPlanar
{        
         struct TwoViewPlanarSDPResult
        {
                // fill by getResult
                double stop_iter; 
                double f_opt; 
                double d_opt; 
                double primal_error;
                double dual_error; 
        
        
                Vector3 dual_point;   // lagrange multipliers
                
                Matrix5 primal_X;  
                
                // modified by getSolutionFromResult
                Vector4 delta_p;   
                
                
                Vector5 eigenvalues_X; 
                
                double f_opt_sol; 
                
        
        };  // end of struct
        
           
        class TwoViewPlanarSDP 
        {
        
        public: 
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                
                TwoViewPlanarSDP(void) {};
                
                TwoViewPlanarSDP(const Matrix3 & H, const Vector3 & p1, const Vector3 & p2, const bool debug = false);  
                
                ~TwoViewPlanarSDP(void) { Problem1_.terminate();  }
                             
                TwoViewPlanarSDPResult getResult(bool debug = false);
                                                      
                void getSolutionFromResult(TwoViewPlanarSDPResult & res, bool debug = false); 
                
                
        private:
        
                       SDPA Problem1_; 
                       int nConstraints_;
                       bool debug_; 
                                                      
        };  // end of class left
               

        
}  // end of namespace PrimalRotPrior
