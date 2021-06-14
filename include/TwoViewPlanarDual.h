#pragma once

#include <eigen3/Eigen/Dense>

// include types
#include "TwoViewPlanarTypes.h"
// include basic fcns
#include "TwoViewPlanarUtils.h"


namespace TwoViewPlanar
{
  
        class TwoViewPlanarClass
        {
                public:
                        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
                        
                        /* Default constructor */
                        TwoViewPlanarClass(void){
                        
                                H_ = Matrix3::Zero(); 
                                U_ = Matrix2::Zero(); 
                                UU_ = Matrix2::Zero(); 
                                h31_ = 1; 
                                h32_ = 1; 
                                s_ = 1; 
                                change_A1_ = 1.0; 
                                
                                max_iter_ = 10;
                                
                                tol_ = 5e-16;   // for constr = 0
                                tol_d_ = 5e-16; // for diff lag mult
                                tol_e_ = 1e-05; 
                
                                time_dec_H_ = 0.0; 
                                };
                
                
                        /* Constructor */
                        TwoViewPlanarClass(Matrix3 & H)
                        {
                                // 1. Decompose E
                                HDecSparseResult res = decomposeHSparse(H); 
                                
                                // 2. Save data
                                H_ = H; 
                                
                                h31_ = res.h31; 
                                h32_ = res.h32; 
                                s_ = res.s;
                                U_ = res.U;
                                UU_ = res.UU;  
                                
                                max_iter_ = 10;
                                
                                change_A1_ = res.change_A1; 
                                
                                tol_ = 5e-16;   // for constr = 0
                                tol_d_ = 5e-16; // for diff lag mult
                                tol_e_ = 1e-05; 
                
                                time_dec_H_ = res.time_dec;
                        }
                        
                        /* Default destructor */
                        ~TwoViewPlanarClass(void){}; 
                        
                        TwoViewPlanarResult correctMatches(const Vector3 & p1, const Vector3 & p2, bool use_affine = false); 
                        
                        /* Print result */ 
                        void printResult(TwoViewPlanarResult & res);
                        
                private:
                        /** Variables **/
                        Matrix3 H_;
                        
                        Matrix2 U_, UU_;
                        
                        double change_A1_; 
                        
                        // elements of T1^TH
                        // 0 0 
                        // a b
                        double h31_; 
                        double h32_; 
                        
                        double s_;   // singular value of Hs2
                        
                        int max_iter_;
                        
                        double tol_;   
        
                        double tol_d_; 
                        
                        double tol_e_; 
                        
                        double time_dec_H_; 
                                                
                        
                        /** Private functions **/                        
                        double _evalGrad12(double l1, double l2, Vector23 & coeff, double norm_bh, double a);

                        double _evalGrad22(double l1, double l2, Vector24 & coeff, double norm_h1_sq, double a);

                        double _evalGrad11(double l1, double l2, Vector24 & coeff, double norm_d1_sq, double a);



                        double _evalConstr2(double l1, double l2, Vector19 & coeff, double c1, double a);

                        double _evalConstr1(double l1, double l2, Vector19 & coeff, double c1, double a);


                        double findZeroNewton(double l_init_1, double l_init_2, 
                                                                          Vector19 & C1_coeff, Vector19 &  C2_coeff, 
                                                                          Vector24 & J11_coeff, Vector24 & J22_coeff, 
                                                                          Vector23 & J12_coeff, 
                                                                          double c1, double c2, 
                                                                          double a, double dot_bh, 
                                                                          double norm_b, double norm_h, 
                                                                          double & l1, double & l2, 
                                                                          double & c1_opt, double & c2_opt);
                                                                          
                                                  
                        
        
        };  // end of main class


}  // end of namespace
