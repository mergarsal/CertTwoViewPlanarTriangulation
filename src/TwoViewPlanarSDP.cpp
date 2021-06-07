#include "TwoViewPlanarSDP.h"


#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

namespace TwoViewPlanar
{



TwoViewPlanarSDP::TwoViewPlanarSDP(const Matrix3& H, const Vector3 & p1, const Vector3 & p2, const bool debug)
{ 
 
                Matrix3 T1 = Matrix3::Zero(), T2 = Matrix3::Zero();
                T1 << 0, 0, 0, 0, 0, -1, 0, 1, 0; 
                T2 << 0, 0, 1, 0, 0, 0, -1, 0, 0; 
                Matrix23 Pro; 
                Pro << 1, 0, 0, 0, 1, 0; 
                
                // we normalize 
                // Matrix3 H = H_original / H_original.trace();
                
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
                               
                
                
     
  Problem1_ = SDPA(); 
  
 
  nConstraints_ = 2 + 1; // 2: minimal + 1: homo.
  int nBlock = 1;    // number blocks
  
  
  // display info from solver
  if (debug) Problem1_.setDisplay(stdout);
  
  
  // All parameteres are renewed
  Problem1_.setParameterType(SDPA::PARAMETER_DEFAULT);
  // Problem1_.printParameters(stdout);

  if (debug) std::cout << "[SDP-PLANAR] Create type solution\n";
  
  Problem1_.inputConstraintNumber(nConstraints_);
  Problem1_.inputBlockNumber(nBlock);
  // bLOCKsTRUCT :: 9(SDP)  6(SDP) 
  Problem1_.inputBlockSize(1,5);
  Problem1_.inputBlockType(1,SDPA::SDP);

  Problem1_.initializeUpperTriangleSpace();
  
  if (debug) std::cout << "[SDP-PLANAR] Define vector b\n";
  //cVECT = {1,0, ..., 0}
  Problem1_.inputCVec(1,1);
  
  // NOTE: needed??
  for (int i=2; i <= nConstraints_; i++) Problem1_.inputCVec(i,0); 


  // ---------  Input F_0  (cost) --------------------

  if (debug) std::cout << "[SDP-PLANAR] Define data matrix\n";
  // 1st block: C
  for (int i = 1; i <= 4; i++)
  {
      // NOTE: SDPA consider the problem as MAX
      Problem1_.inputElement(0, 1, i, i, - 1 );    
  }

  

  /* Define constraints */
  
  // auxiliary variables
  size_t d11=1, d12=2, d21=3, d22=4;
  int y = 5; 
  
  /* Homogeneous constraint */ 
  // # constraint, block, i1, j1, value 
  Problem1_.inputElement(1, 1, y, y, 1);
  
  
  if (debug) std::cout << "[SDP-PLANAR] Define constraint c1\n";
  // ---------  Input F_2  --------------------

  
  size_t id_const = 2; 
  // # constraint, block, i1, j1, value 
  /* First set */ 
  Problem1_.inputElement(id_const, 1, d11, d21, Hs1(0,0));
  Problem1_.inputElement(id_const, 1, d12, d21, Hs1(0,1));
  Problem1_.inputElement(id_const, 1, d22, d11, Hs1(1,0));
  Problem1_.inputElement(id_const, 1, d22, d12, Hs1(1,1));
  Problem1_.inputElement(id_const, 1, y, d11, k1(0));
  Problem1_.inputElement(id_const, 1, y, d12, k1(1));
  Problem1_.inputElement(id_const, 1, y, d21, f1(0));
  Problem1_.inputElement(id_const, 1, y, d22, f1(1));
  Problem1_.inputElement(id_const, 1, y, y, c1);

  id_const++;
  Problem1_.inputElement(id_const, 1, d11, d21, Hs2(0,0));
  Problem1_.inputElement(id_const, 1, d12, d21, Hs2(0,1));
  Problem1_.inputElement(id_const, 1, d22, d11, Hs2(1,0));
  Problem1_.inputElement(id_const, 1, d22, d12, Hs2(1,1));
  Problem1_.inputElement(id_const, 1, y, d11, k2(0));
  Problem1_.inputElement(id_const, 1, y, d12, k2(1));
  Problem1_.inputElement(id_const, 1, y, d21, f2(0));
  Problem1_.inputElement(id_const, 1, y, d22, f2(1));
  Problem1_.inputElement(id_const, 1, y, y, c2);

  
    
}



TwoViewPlanarSDPResult TwoViewPlanarSDP::getResult(bool debug){


  /* Solve problem now */ 
  if (debug) std::cout << "[SDP-PLANAR] Getting result\n";
  Problem1_.initializeUpperTriangle();  // early termination
  if (debug) std::cout << "[SDP-PLANAR] Init solver\n";
  Problem1_.initializeSolve();
  
  // if necessary, dump input data and initial point
  if (debug) std::cout << "[SDP-PLANAR] Solve!\n";
  Problem1_.solve();
  
  //
  if (debug) std::cout << "[SDP-PLANAR] Saving results\n";
  TwoViewPlanarSDPResult res = TwoViewPlanarSDPResult(); 
  res.stop_iter = Problem1_.getIteration(); 
  res.d_opt = - Problem1_.getPrimalObj(); 
  res.f_opt = - Problem1_.getDualObj();
  res.primal_error = Problem1_.getPrimalError(); 
  res.dual_error = Problem1_.getDualError(); 
  
  
  
 
  
  if (debug) std::cout << "[SDP-PLANAR] Saving lagrange multipliers\n";
  auto dual_point = Problem1_.getResultXVec(); 
  for (int i=0; i < nConstraints_; i++)
        res.dual_point(i) = dual_point[i]; 
  
  

  if (debug) std::cout << "[SDP-PLANAR] Saving blocks for solutions\n";
  auto X = Problem1_.getResultYMat(1);  

  // std::cout << "Printing Xe\n"; 
    
  for (int i=0; i<5; ++i) {
    for (int j=0; j<5; ++j) {
      res.primal_X(i, j) = X[i * 5 + j];  // 5 x 5
      }
      }
      
      
   
  return res; 
}


/** 
Get results for this problem
**/  
  

void TwoViewPlanarSDP::getSolutionFromResult(TwoViewPlanarSDPResult & res, bool debug){



    auto max_index = [] (const Eigen::Matrix<double, Eigen::Dynamic, 1> & e_vector, 
                                                        double val_rot) -> int
                {
                int i = 0;  
                while (i < e_vector.rows())
                {
                        if (e_vector[i] == val_rot) 
                                break;   
                        else i++; 
                }
                return i;
                };
                
                
                
                
                
        // Extract solution from blocks
        // This function modifies some of the fields in res
        // E_opt, t_opt, q_opt
        Matrix5 X = res.primal_X; 
        
       if (debug)  std::cout << "Extracting data from SDP solution\n"; 
         
         
        // Extract essential matrix
        
        Eigen::SelfAdjointEigenSolver<Matrix5> eigen_solver_M(X);
                
        Vector5 eigenvalues_X = eigen_solver_M.eigenvalues().real();
        
        double val_rot = eigenvalues_X.maxCoeff();  
                

                int mm = max_index(eigenvalues_X, val_rot); 
                Vector5 sol_x = eigen_solver_M.eigenvectors().col(mm);
                
                if (debug)
                {
                        std::cout << "Eigenvalues of Xe:\n" << eigenvalues_X.transpose(); 
                        std::cout << "\nMaximum eigenvalue = " << val_rot; 
                        std::cout << " with index " << mm << std::endl;
                        std::cout << "Eigenvector:\n" << sol_x.transpose() << std::endl;
                
                }
                
                // reshape ee 
                Vector4 delta_p;  
                // fill up expanded e
                delta_p << sol_x(0) / sol_x(4), sol_x(1) / sol_x(4), sol_x(2) / sol_x(4), sol_x(3) / sol_x(4);
                
                res.delta_p = delta_p; 
                res.eigenvalues_X = eigenvalues_X; 
                res.f_opt_sol = delta_p.dot(delta_p);                

        return;
}







}  // end namesapce 
