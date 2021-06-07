#include "triangulationPlanar.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>



#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <Eigen/QR>  

#include <vector>

int rpoly(double *op, int degree, double *zeror, double *zeroi, int info[] );


namespace TwoViewPlanar
{


        void computeHomography(Vector3 & n, double d_plane, 
                               Matrix3 & R1, Matrix3 & R2, 
                               Vector3 & t1, Vector3 & t2, 
                               double focal_length, double size_img, 
                               Matrix3 & H, Matrix3 & Hp)
                               
                               
        {
                Matrix3 Rrel = R2 * R1.transpose(); 
                Vector3 Trel = t2 - Rrel * t1; 
                
                Vector3 m = R1 * n; 
                double d1 = d_plane + m.dot(t1);
                
                H.setZero(); 
                Hp.setZero(); 
                H = Rrel + Trel * m.transpose() / d1; 
                
                Matrix3 K = Matrix3::Zero(); 
                K(0,0) = focal_length; 
                K(1,1) = focal_length; 
                K(0,2) = size_img / 2; 
                K(1,2) = size_img / 2; 
                K(2,2) = 1.0; 
                Hp = K * H * K.inverse();     
        
        
        }

        
        
        int correctMatchesPoly (const Matrix3& H,
                                const Vector2& point1,
                                const Vector2& point2,
                                Vector2* corrected_point1,
                                Vector2* corrected_point2)
                            
                            {
                                   
                            
                                    /* 1. Create transformations 1: translations */     
                                    
                                 
                                    Matrix3 Li = Matrix3::Identity(), Lp = Matrix3::Identity(); 
                                    
                                    // NOTE: here L is the inverse, so +
                                    Li(0, 2) = point1(0); 
                                    Li(1, 2) = point1(1);
                                    
                                    Lp(0, 2) = -point2(0); 
                                    Lp(1, 2) = -point2(1);
                                    
                                   
                                    
                                    // 2. Create B = Lp * H * Li
                                    Matrix3 B = Lp*H*Li;
                                    
                                    
                                    // 4. Compute Q = B * Ri = B * R^T
                                    // alpha = arctan(-h8/h7)
                                    double sr, cr, norm_h87;
                                    norm_h87 = std::sqrt(H(2, 0)*H(2, 0) + H(2, 1)*H(2, 1)); 
                                    Matrix3 R = Matrix3::Identity(); 
                                    if (norm_h87 > 1e-10)
                                    {
                                            sr = H(2, 1) / norm_h87;
                                            cr = H(2, 0) / norm_h87; 
                                            R(0, 0) = cr;
                                            R(0, 1) = sr; 
                                            R(1, 1) = cr; 
                                            R(1, 0) = -sr; 
                                    }                                 
                                    
                                    
                                    
                                    Matrix3 Q = B * R.transpose(); 
                                    
                                   
                                    
                                    // Apply second rotation 
                                    // Let Rp so Qp = Rp * Q 
                                    Matrix3 Rp = Matrix3::Identity(); 
                                    double srp, crp, norm_q14; 
                                    norm_q14 = std::sqrt(Q(0,0)*Q(0,0) + Q(1,0) * Q(1, 0));
                                    if (norm_q14 > 1e-10)
                                    {
                                            srp = Q(1, 0) / norm_q14; 
                                            crp = Q(0, 0) / norm_q14; 
                                            Rp(0, 0) = crp;
                                            Rp(0, 1) = srp; 
                                            Rp(1, 1) = crp; 
                                            Rp(1, 0) = -srp; 
                                    }                                 
                                    
                                    
                                    
                                    Matrix3 Qp = Rp * Q; 
                                    
                                  
                                    
                                    
                                    // Get some values
                                    double t, r;
                                    double q1, q2, q3, q5, q6, q7, q9; 
                                    q1 = Qp(0,0); 
                                    q2 = Qp(0,1); 
                                    q3 = Qp(0,2); 
                                    q5 = Qp(1,1); 
                                    q6 = Qp(1,2);
                                    q7 = Qp(2,0); 
                                    q9 = Qp(2,2); 
                                    
                                    t = q3 * q5 - q2 * q6; 
                                    r = q2*q2 + q5*q5 + q9 * q9;   
                                    
                                   
                                   // 6. Form the polynomial 
                                   // g(t) = k0 * t^8 + k1 * t^7 + k2*t^6 + k3*t^5 + k4*t^4 + k5*t^3 + k6*t^2 + k7*t + k8
                                   double k0, k1, k2, k3, k4, k5, k6, k7, k8; 
                                   
                                   k8 = q9*q9*q9 * ((-q3*q3 -q6*q6)*q7*q9+q1*q3*r) + q1*q5*q9*r*t -q7*(q9*q9 +r)*t*t; 
                                   
                                   k7 = -4*q3*q3*q7*q7*q9*q9*q9 - 4*q6*q6*q7*q7*q9*q9*q9 + 3*q1*q3*q7*q9*q9*r + q9*r*(q1*q1*(q5*q5+q9*q9)+q9*q9*r) - q1*q5*q7*r*t-4*q7*q7*q9*t*t; 
                                   
                                   k6 = q7*(q9*(-6*(q3*q3+q6*q6)*q7*q7*q9 - q1*q3*q7*(q9*q9-3*r)+4*q9*q9*q9*r +3*q9*r*r + q1*q1*q9*(q5*q5+q9*q9+3*r))-5*q1*q5*q7*q9*t-2*q7*q7*t*t);
                                   
                                   k5 = q7*q7*(q9*(-4*(q3*q3+q6*q6)*q7*q7+4*q9*q9*q9*q9+14*q9*q9*r+3*r*r)+q1*q1*(-(q5*q5*q9)+3*q9*(q9*q9+r))+q1*q7*(q3*(-3*q9*q9+r)-3*q5*t));
                                   
                                   k4 = q7*q7*q7*((-q3*q3-q6*q6)*q7*q7-3*q1*q3*q7*q9+16*q9*q9*q9*q9+18*q9*q9*r+r*r+q1*q1*(-q5*q5+3*q9*q9+r)); 
                                   
                                   k3 = q7*q7*q7*q7*(-(q1*q3*q7)+q1*q1*q9+25*q9*q9*q9+10*q9*r); 
                                   
                                   k2 = q7*q7*q7*q7*q7*(19*q9*q9+2*r); 
                                   
                                   k1 = 7*q7*q7*q7*q7*q7*q7*q9; 
                                   
                                   k0 = q7*q7*q7*q7*q7*q7*q7;
                                   
                                   
                                   int n_coeff = 8; 
                                   
                                   if (std::abs(k0) < 1e-10)
                                   {
                                        if (std::abs(k1) > 1e-10)       n_coeff = 7;
                                        else
                                        {
                                                if (std::abs(k2) > 1e-10)       n_coeff = 6;
                                                else
                                                {
                                                        if (std::abs(k3) > 1e-10)     n_coeff = 5;
                                                        else
                                                        {
                                                                if (std::abs(k4) > 1e-10)        n_coeff = 4;
                                                                else
                                                                {
                                                                        if (std::abs(k5) > 1e-10)    n_coeff = 3;
                                                                        else
                                                                        {
                                                                                if (std::abs(k6) > 1e-10)
                                                                                        n_coeff = 2;
                                                                                else
                                                                                        if (std::abs(k7) > 1e-10)
                                                                                                n_coeff = 1;
                                                                                        else
                                                                                                return -1; 
                                                                        }
                                                                }
                                                        }
                                                }                                        
                                        }
                                   }
                                   // 7. Call solver 
                                   double * poly_coeff = new double [n_coeff+1];
                                   
                                   int id_c = n_coeff;
                                   poly_coeff[id_c] = k8; 
                                   poly_coeff[--id_c] = k7;  
                                   if (n_coeff > 1)      
                                   {                          
                                           poly_coeff[--id_c] = k6;
                                           if (n_coeff > 2)
                                           {
                                                   poly_coeff[--id_c] = k5;
                                                   if (n_coeff > 3)
                                                   {
                                                           poly_coeff[--id_c] = k4;
                                                           if (n_coeff > 4)
                                                           {
                                                                   poly_coeff[--id_c] = k3;
                                                                   if (n_coeff > 5)
                                                                   {
                                                                           poly_coeff[--id_c] = k2;
                                                                           if (n_coeff > 6)
                                                                           {
                                                                                   poly_coeff[--id_c] = k1; 
                                                                                   if (n_coeff > 7)  poly_coeff[--id_c] = k0; 
                                                                           }
                                                                    }
                                                           }
                                                    }
                                            }
                                    }
                                   
                                     
                                         
                                         
                                   double * r_zeros = new double [n_coeff], * i_zeros = new double [n_coeff];
                                   int * info = new int [n_coeff+1];
                                   
                                   
                                   
                               
                                   int n_roots = rpoly(poly_coeff, n_coeff, r_zeros, i_zeros, info);
                                   
                                   
                                   
                                  
                                   
                                   // 8. Select optimal solution 
                                   // Evaluate the cost function s(t) at the real part of the 6 roots
                                   double s_val = 1000000000000; 
                                   double xh_min, xhp_min; 
                                   double yh_min, yhp_min; 
                                   
                                   for (int i=0; i < n_roots; i++)
                                   {
                                        double xh = r_zeros[i]; 
                                        double yh = - (q2*q3+q5*q6 + q1*q2*xh) / (q2*q2+q5*q5+q9*q9+2*q7*q9*xh+q7*q7*xh*xh);
                                        
                                        double xhp = (q1*xh + q2*yh + q3) / (q7*xh + q9) ;   
                                        double yhp = (q5*yh+q6) / (q7*xh+q9);
                                        double s = xh*xh+yh*yh + xhp*xhp + yhp*yhp; 
                                         
                                         
                                        
                                         if (s < s_val) {
                                                s_val = s;
                                                xh_min = xh;
                                                yh_min = yh; 
                                                yhp_min = yhp; 
                                                xhp_min = xhp;
                                         }
                                   
                                   }
                                                                      
                                   
                                   
                                  
                                   
                                   
                                   // 9. Recover points
                                   Vector3 temp = Vector3::Zero(); 
                                   temp(0) = xh_min;
                                   temp(1) = yh_min; 
                                   temp(2) = 1; 
                                                                      
                                   
                 
                                   Vector3 temp_rt = Li*R.transpose()*temp; 
                                   temp_rt /= temp_rt(2); 
                                   
                                   
                                   
                                   *corrected_point1 << temp_rt(0),  temp_rt(1); 
                                   
                                   
                                   
                                   
                                   Vector3 temp_p = Vector3::Zero(); 
                                   temp_p(0) = xhp_min; 
                                   temp_p(1) = yhp_min; 
                                   temp_p(2) = 1;                                    
                                   
                                   Vector3 temp_p_rot = Lp.inverse() * Rp.transpose() * temp_p; 
                                   temp_p_rot /= temp_p_rot(2); 
                                                                      
                                   *corrected_point2 << temp_p_rot(0), temp_p_rot(1); 
                               
                                   
                                   return n_roots;   
                            }  // end of correctMatchesPoly

        
        
        
        // returns number of iterations
        int correctMatchesKanatani (const Matrix3& H,
                                const Vector2& point1,
                                const Vector2& point2,
                                const double f0, 
                                Vector2* corrected_point1,
                                Vector2* corrected_point2)
                                {
                                
                                double E = 1000, E0 = 10000;  // initial (current) error 
                                
                                Eigen::Matrix<double, 4, 1> p, ph, ptilde; 
                                p << point1, point2; 
                                ph = p; 
                                ptilde.setZero(); 
                                
                                
                                Eigen::Matrix<double, 9, 4> T1, T2, T3; 
                                Eigen::Matrix<double, 9, 1> e1, e2, e3; 
                                
                                Eigen::Matrix<double, 9, 9> V11, V12, V13, V22, V23, V33, V21, V31, V32; 
                                
                                Eigen::Matrix<double, 3, 3> W, Dw, iW;  // matrix W and its inverse 
                                Vector3 eigen_W; 
                                
                                Matrix3 Ht = H.transpose(); 
                                Eigen::Matrix<double, 9, 1> h = Eigen::Map<const Eigen::Matrix<double, 9, 1>>(Ht.data(), H.size());                               
                                
                                
                           
                                
                                bool conv = false;  // flag for convergence
                                
                                int iter = 0;                                 
                                for (iter = 0; iter < 500; iter++)
                                {
                                        // 2. Fill T's
                                        T1.setZero(); 
                                        T1(3, 0) = - f0; 
                                        T1(4, 1) = - f0; 
                                        T1(6, 3) = ph(0); 
                                        T1(6, 0) = ph(3); 
                                        T1(7, 1) = ph(3); 
                                        T1(7, 3) = ph(1); 
                                        T1(8, 3) = f0; 
                                        
                                        T2.setZero(); 
                                        T2(0, 0) = f0; 
                                        T2(1, 1) = f0; 
                                        T2(6, 0) = -ph(2); 
                                        T2(6, 2) = -ph(0); 
                                        T2(7, 1) = -ph(2); 
                                        T2(7, 2) = -ph(1); 
                                        T2(8, 3) = -f0;
                                        
                                        T3.setZero();
                                        T3(0, 0) = -ph(3); 
                                        T3(0, 3) = -ph(0); 
                                        T3(1, 1) = -ph(3); 
                                        T3(1, 3) = -ph(1); 
                                        T3(2, 3) = -f0; 
                                        T3(3, 0) = ph(2); 
                                        T3(3, 2) = ph(0); 
                                        T3(4, 1) = ph(2); 
                                        T3(4, 2) = ph(1);
                                        T3(5, 2) = f0; 
                                                                               
                                        double x, y, xp, yp; 
                                        x = ph(0); 
                                        y = ph(1); 
                                        xp = ph(2); 
                                        yp = ph(3); 
                                        // 3. Compute E's
                                        e1 << 0, 0, 0, -f0*x, -f0*y, -f0*f0, x*yp, y*yp, f0*yp; 
                                        e1 += T1 * ptilde; 
                                        
                                        e2 << f0*x, f0*y, f0*f0, 0, 0, 0, -x*xp, -y*xp, -f0*xp;
                                        e2 += T2 * ptilde; 
                                        
                                        e3 << -x*yp, -y*yp, -f0*yp, x*xp, y*xp, f0*xp, 0, 0, 0; 
                                        e3 += T3 * ptilde; 
                                        
                                        
                                        // 4. Compute the 9x9 matrices Vo(k, l)
                                        V11 = T1 * T1.transpose(); 
                                        V12 = T1 * T2.transpose(); 
                                        V13 = T1 * T3.transpose(); 
                                        V22 = T2 * T2.transpose(); 
                                        V23 = T2 * T3.transpose(); 
                                        V33 = T3 * T3.transpose(); 
                                        V21 = T2 * T1.transpose(); 
                                        V31 = T3 * T1.transpose(); 
                                        V32 = T3 * T2.transpose(); 
                                        
                                        // 5. Compute the 3x3 matrix W
                                        W.setZero(); 
                                        W(0, 0) = h.dot(V11 * h); 
                                        W(0, 1) = h.dot(V12 * h); 
                                        W(0, 2) = h.dot(V13 * h); 
                                        W(1, 0) = h.dot(V21 * h); 
                                        W(1, 1) = h.dot(V22 * h); 
                                        W(1, 2) = h.dot(V23 * h); 
                                        W(2, 0) = h.dot(V31 * h); 
                                        W(2, 1) = h.dot(V32 * h); 
                                        W(2, 2) = h.dot(V33 * h); 
                                        

                                        
                                        // We are suppoosed to compute the inverse with rank 2
                                        // (the smallest eigenvalue is replaced by 0 in its spectral decomposition).
                                      
                                        
                                        
                                        Eigen::SelfAdjointEigenSolver<Matrix3> eigensolver(W);
                                        
                                        
                                        eigen_W = eigensolver.eigenvalues(); 
                                        Dw.setZero(); 
                                        Dw(1,1) = eigen_W(1); 
                                        Dw(2,2) = eigen_W(2);
                                        
                                       
                                        
                                        // Update
                                        W = eigensolver.eigenvectors() * Dw * eigensolver.eigenvectors().transpose(); 
                                        
                                     
                                        // Option 1: 
                                        iW = W.completeOrthogonalDecomposition().pseudoInverse();
                                 

                                        // 6. Update ptilde as phat as follows 
                                        // k=1, l= 1
                                        ptilde =  iW(0, 0) * e1.dot(h) * T1.transpose() * h;
                                        // k = 1, l = 2
                                        ptilde += iW(0, 1) * e2.dot(h) * T1.transpose() * h; 
                                        // k = 1, l = 3
                                        ptilde += iW(0, 2) * e3.dot(h) * T1.transpose() * h;
                                        
                                        // k = 2, l = 1
                                        ptilde += iW(1, 0) * e1.dot(h) * T2.transpose() * h;
                                        // k = 2, l = 2
                                        ptilde += iW(1, 1) * e2.dot(h) * T2.transpose() * h;
                                        // k = 2, l = 3
                                        ptilde += iW(1, 2) * e3.dot(h) * T2.transpose() * h;
                                        
                                        // k = 3, l = 1
                                        ptilde += iW(2, 0) * e1.dot(h) * T3.transpose() * h;
                                        // k = 3, l = 2
                                        ptilde += iW(2, 1) * e2.dot(h) * T3.transpose() * h;
                                        // k = 3, l = 3
                                        ptilde += iW(2, 2) * e3.dot(h) * T3.transpose() * h;
                                        
                                        // update ph
                                        ph = p - ptilde;
                                        
                                        // 7. Evaluate the reprojection error E = norm(ptild)^2
                                        // a. If E approx Eo, then return ph ans p_corrected and stop
                                        E = ptilde.dot(ptilde); 
                                       
                                        // else, E0 <- E and go back to iter++
                                        if (std::abs(E0 - E) < 5e-14)
                                        {
                                                conv = true;  
                                                break;                                      
                                        }
                                        
                                        
                                        // else
                                        E0 = E;
                                        // go to next iter
                                        
                                
                                }  // end of loop
                                
                                if (conv)
                                {
                                        *corrected_point1 << ph(0), ph(1); 
                                        *corrected_point2 << ph(2), ph(3);  
                                
                                }
                                
                                return ++iter; 
                                
                                
                                }  // end of correctMatchesKanatani
                                
                                
        


                           int FindObsHS (const Eigen::Matrix<double, 3, 3>& ematrix,
                            const Eigen::Matrix<double, 2, 1>& point1,
                            const Eigen::Matrix<double, 2, 1>& point2,
                            Eigen::Matrix<double, 2, 1>* corrected_point1,
                            Eigen::Matrix<double, 2, 1>* corrected_point2)
                            
                            {
                                   
                            
                                    /* 1. Create transformations 1: translations */                                    
                                    // NOTE: these are the inverted of the transformations
                                    // that's why the sign are +
                                    
                                    
                                    Eigen::Matrix<double, 3, 3> T1i = Eigen::Matrix<double, 3, 3>::Identity(), T2i = Eigen::Matrix<double, 3, 3>::Identity(); 
                                    T1i(0, 2) = point1(0); 
                                    T1i(1, 2) = point1(1);
                                    
                                    T2i(0, 2) = point2(0); 
                                    T2i(1, 2) = point2(1);
                                    
                                    
                                    // 2. Create T2i*F*T1i
                                    Eigen::Matrix<double, 3, 3> TET = T2i.transpose() * ematrix.transpose() * T1i; 
                                    
                                    // 3. Compute epipoles                                     
                                    Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd(TET, Eigen::ComputeFullU | Eigen::ComputeFullV);
                                    
                                    // a. Left epipole
                                    Eigen::Matrix<double, 3, 1> e2 = svd.matrixU().col(2);
                                    // compute scale. we do not assume that the epipole is within the image
                                    double scale; 
                                    scale = sqrt(e2(0)*e2(0) + e2(1) * e2(1));
                                    e2 /= scale;         
                                        
                                    // b. right epipole
                                    Eigen::Matrix<double, 3, 1> e1 = svd.matrixV().col(2);
                                    scale = sqrt(e1(0)*e1(0) + e1(1) * e1(1));
                                    e1 /= scale;  
                                    
                                    // 4. Compute RTETR
                                    Eigen::Matrix<double, 3, 3> R1 = Eigen::Matrix<double, 3, 3>::Identity(), R2 = Eigen::Matrix<double, 3, 3>::Identity(); 
                                    R1(0, 0) = e1(0); 
                                    R1(0, 1) = e1(1); 
                                    R1(1, 0) = -e1(1); 
                                    R1(1, 1) = e1(0); 
                                    
                                    R2(0, 0) = e2(0); 
                                    R2(0, 1) = e2(1); 
                                    R2(1, 0) = -e2(1); 
                                    R2(1, 1) = e2(0); 
                            
                                    
                                    Eigen::Matrix<double, 3, 3> RER = R2 * TET * R1.transpose(); 
                                    
                                    
                                    // 5. Save some vars
                                    double f1, f2, a, b, c, d;
                                    f1 = e1(2); 
                                    f2 = e2(2); 
                                    a = RER(1, 1); 
                                    b = RER(1, 2);
                                    c = RER(2, 1); 
                                    d = RER(2, 2);  
                                    
                                   
                                   // 6. Form the polynomial g(t) = k6*t^6 + k5*t^5 + k4*t^4 + k3*t^3 + k2*t^2 + k1*t + k0
                                   double k0, k1, k2, k3, k4, k5, k6; 
                                   
                                   k0 = -a*d*d*b+b*b*c*d;
                                   k1 = +f2*f2*f2*f2*d*d*d*d+b*b*b*b+2*b*b*f2*f2*d*d-a*a*d*d+b*b*c*c ;
                                   k2 = +4*a*b*b*b+4*b*b*f2*f2*c*d+4*f2*f2*f2*f2*c*d*d*d-a*a*d*c+b*c*c*a+4*a*b*f2*f2*d*d-2*a*d*d*f1*f1*b+2*b*b*c*f1*f1*d; 
                                   k3 = +6*a*a*b*b+6*f2*f2*f2*f2*c*c*d*d+2*b*b*f2*f2*c*c+2*a*a*f2*f2*d*d-2*a*a*d*d*f1*f1+2*b*b*c*c*f1*f1+8*a*b*f2*f2*c*d; 
                                   k4 = +4*a*a*a*b+2*b*c*c*f1*f1*a+4*f2*f2*f2*f2*c*c*c*d+4*a*b*f2*f2*c*c+4*a*a*f2*f2*c*d-2*a*a*d*f1*f1*c-a*d*d*f1*f1*f1*f1*b+b*b*c*f1*f1*f1*f1*d; 
                                   k5 = +f2*f2*f2*f2*c*c*c*c+2*a*a*f2*f2*c*c-a*a*d*d*f1*f1*f1*f1+b*b*c*c*f1*f1*f1*f1+a*a*a*a; 
                                   k6 = +b*c*c*f1*f1*f1*f1*a-a*a*d*f1*f1*f1*f1*c; 
                                   
                                   // 7. Call solver 
                                   double * poly_coeff = new double [7];
                                   poly_coeff[6] = k0; 
                                   poly_coeff[5] = k1;
                                   poly_coeff[4] = k2;
                                   poly_coeff[3] = k3;
                                   poly_coeff[2] = k4;
                                   poly_coeff[1] = k5;
                                   poly_coeff[0] = k6;
                                   
                                   double * r_zeros = new double [6], * i_zeros = new double [6];
                                   int * info = new int [7];
                                   
                                   int n_roots = rpoly(poly_coeff, 6, r_zeros, i_zeros, info);
                                   
                                   
                                   double s_val = 1./(f1*f1) + (c*c)/(a*a+f2*f2*c*c);
                                   double t_min; 
                                   for (int i=0; i < n_roots; i++)
                                   {
                                        double t = r_zeros[i]; 
                                        double s = (t*t)/(1 + f1*f1*t*t) + ((c*t + d)*(c*t + d))/((a*t + b)*(a*t + b) + f2*f2*(c*t + d)*(c*t + d));
                                        
                                         if (s < s_val) {
                                                s_val = s;
                                                t_min = t;
                                         }
                                   
                                   }
                                   
                                   // 9. Recover points
                                   Eigen::Matrix<double, 3, 1> temp = Eigen::Matrix<double, 3, 1>::Zero(); 
                                   temp(0) = t_min*t_min*f1; 
                                   temp(1) = t_min; 
                                   temp(2) = t_min*t_min*f1*f1+1;  
                                   Eigen::Matrix<double, 3, 1> temp_rt = T1i*R1.transpose()*temp; 
                                   temp_rt /= temp_rt(2); 
                                   *corrected_point1 << temp_rt(0),  temp_rt(1); 
                                   
                                   temp.setZero(); 
                                   temp_rt.setZero(); 
                                   temp(0) = f2*(c*t_min+d)*(c*t_min+d);
                                   temp(1) = -(a*t_min+b)*(c*t_min+d);
                                   temp(2) = f2*f2*(c*t_min+d)*(c*t_min+d) + (a*t_min+b)*(a*t_min+b);
                                   
                                   temp_rt = T2i*R2.transpose()*temp; 
                                   temp_rt /= temp_rt(2); 
                                   *corrected_point2 << temp_rt(0), temp_rt(1); 
                                   
                                   return n_roots;   
                            }




}  // end of namespace










/** From: http://www.crbond.com/download/misc/rpoly.cpp **/
/*      rpoly.cpp -- Jenkins-Traub real polynomial root finder.
 *
 *      (C) 2002, C. Bond.  All rights reserved.
 *
 *      Translation of TOMS493 from FORTRAN to C. This
 *      implementation of Jenkins-Traub partially adapts
 *      the original code to a C environment by restruction
 *      many of the 'goto' controls to better fit a block
 *      structured form. It also eliminates the global memory
 *      allocation in favor of local, dynamic memory management.
 *
 *      The calling conventions are slightly modified to return
 *      the number of roots found as the function value.
 *
 *      INPUT:
 *      op - double precision vector of coefficients in order of
 *              decreasing powers.
 *      degree - integer degree of polynomial
 *
 *      OUTPUT:
 *      zeror,zeroi - output double precision vectors of the
 *              real and imaginary parts of the zeros.
 *
 *      RETURN:
 *      returnval:   -1 if leading coefficient is zero, otherwise
 *                  number of roots found. 
 */


void quad(double a,double b1,double c,double *sr,double *si,
        double *lr,double *li);
void fxshfr(int l2, int *nz);
void quadit(double *uu,double *vv,int *nz);
void realit(double sss, int *nz, int *iflag);
void calcsc(int *type);
void nextk(int *type);
void newest(int type,double *uu,double *vv);
void quadsd(int n,double *u,double *v,double *p,double *q,
        double *a,double *b);
double *p,*qp,*k,*qk,*svk;
double sr,si,u,v,a,b,c,d,a1,a2;
double a3,a6,a7,e,f,g,h,szr,szi,lzr,lzi;
double eta,are,mre;
int n,nn,nmi,zerok;
static int itercnt;

int rpoly(double *op, int degree, double *zeror, double *zeroi, int info[] ) 
{
    double t,aa,bb,cc,*temp,factor,rot;
    double *pt;
    double lo,max,min,xx,yy,cosr,sinr,xxx,x,sc,bnd;
    double xm,ff,df,dx,infin,smalno,base;
    int cnt,nz,i,j,jj,l,nm1,zerok;
    long sec;

    
    sec = clock();

/*  The following statements set machine constants. */
    base = 2.0;
    eta = 2.22e-16;
    infin = 3.4e38;
    smalno = 1.2e-38;

    are = eta;
    mre = eta;
    lo = smalno/eta;
/*  Initialization of constants for shift rotation. */        
    xx = sqrt(0.5);
    yy = -xx;
    rot = 94.0;
    rot *= 0.017453293;
    cosr = cos(rot);
    sinr = sin(rot);
    n = degree;
/*  Algorithm fails of the leading coefficient is zero. */
    if (op[0] == 0.0) return -1;
/*  Remove the zeros at the origin, if any. */
    while (op[n] == 0.0) {
        j = degree - n;
        zeror[j] = 0.0;
        zeroi[j] = 0.0;
        n--;
    }
    if (n < 1) return degree;
/*
 *  Allocate memory here
 */
    temp = new double [degree+1];
    pt = new double [degree+1];
    p = new double [degree+1];
    qp = new double [degree+1];
    k = new double [degree+1];
    qk = new double [degree+1];
    svk = new double [degree+1];
/*  Make a copy of the coefficients. */
    for (i=0;i<=n;i++)
        p[i] = op[i];
/*  Start the algorithm for one zero. */
_40:        
    itercnt = 0;
    if (n == 1) {
        zeror[degree-1] = -p[1]/p[0];
        zeroi[degree-1] = 0.0;
        n -= 1;
        if( info != NULL )
           info[ degree ] = 0;

        goto _99;
    }
/*  Calculate the final zero or pair of zeros. */
    if (n == 2) {
        quad(p[0],p[1],p[2],&zeror[degree-2],&zeroi[degree-2],
            &zeror[degree-1],&zeroi[degree-1]);
        n -= 2;
        if( info != NULL )
           info[ degree ] = info[ degree - 1] = 0;
        goto _99;
    }
/*  Find largest and smallest moduli of coefficients. */
    max = 0.0;
    min = infin;
    for (i=0;i<=n;i++) {
        x = fabs(p[i]);
        if (x > max) max = x;
        if (x != 0.0 && x < min) min = x;
    }
/*  Scale if there are large or very small coefficients.
 *  Computes a scale factor to multiply the coefficients of the
 *  polynomial. The scaling si done to avoid overflow and to
 *  avoid undetected underflow interfering with the convergence
 *  criterion. The factor is a power of the base.
 */
    sc = lo/min;
    if (sc > 1.0 && infin/sc < max) goto _110;
    if (sc <= 1.0) {
        if (max < 10.0) goto _110;
        if (sc == 0.0)
            sc = smalno;
    }
    l = (int)(log(sc)/log(base) + 0.5);
    factor = pow(base*1.0,l);
    if (factor != 1.0) {
        for (i=0;i<=n;i++) 
            p[i] = factor*p[i];     /* Scale polynomial. */
    }
_110:
/*  Compute lower bound on moduli of roots. */
    for (i=0;i<=n;i++) {
        pt[i] = (fabs(p[i]));
    }
    pt[n] = - pt[n];
/*  Compute upper estimate of bound. */
    x = exp((log(-pt[n])-log(pt[0])) / (double)n);
/*  If Newton step at the origin is better, use it. */        
    if (pt[n-1] != 0.0) {
        xm = -pt[n]/pt[n-1];
        if (xm < x)  x = xm;
    }
/*  Chop the interval (0,x) until ff <= 0 */
    while (1) {
        xm = x*0.1;
        ff = pt[0];
        for (i=1;i<=n;i++) 
            ff = ff*xm + pt[i];
        if (ff <= 0.0) break;
        x = xm;
    }
    dx = x;
/*  Do Newton interation until x converges to two 
 *  decimal places. 
 */
    while (fabs(dx/x) > 0.005) {
        ff = pt[0];
        df = ff;
        for (i=1;i<n;i++) { 
            ff = ff*x + pt[i];
            df = df*x + ff;
        }
        ff = ff*x + pt[n];
        dx = ff/df;
        x -= dx;
        itercnt++;
    }
    bnd = x;
/*  Compute the derivative as the initial k polynomial
 *  and do 5 steps with no shift.
 */
    nm1 = n - 1;
    for (i=1;i<n;i++)
        k[i] = (double)(n-i)*p[i]/(double)n;
    k[0] = p[0];
    aa = p[n];
    bb = p[n-1];
    zerok = (k[n-1] == 0);
     
    for(jj=0;jj<5;jj++) {
        itercnt++;
        cc = k[n-1];
        if (!zerok) {
/*  Use a scaled form of recurrence if value of k at 0 is nonzero. */             
            t = -aa/cc;
            for (i=0;i<nm1;i++) {
                j = n-i-1;
                k[j] = t*k[j-1]+p[j];
            }
            k[0] = p[0];
            zerok = (fabs(k[n-1]) <= fabs(bb)*eta*10.0);
        }
        else {
/*  Use unscaled form of recurrence. */
            for (i=0;i<nm1;i++) {
                j = n-i-1;
                k[j] = k[j-1];
            }
            k[0] = 0.0;
            zerok = (k[n-1] == 0.0);
        }
    } 
/*  Save k for restarts with new shifts. */
    for (i=0;i<n;i++) 
        temp[i] = k[i];
/*  Loop to select the quadratic corresponding to each new shift. */ 
    for (cnt = 0;cnt < 20;cnt++) {
/*  Quadratic corresponds to a double shift to a            
 *  non-real point and its complex conjugate. The point
 *  has modulus bnd and amplitude rotated by 94 degrees
 *  from the previous shift.
 */ 
        xxx = cosr*xx - sinr*yy;
        yy = sinr*xx + cosr*yy;
        xx = xxx;
        sr = bnd*xx;
        si = bnd*yy;
        u = -2.0 * sr;
        v = bnd;
        fxshfr(20*(cnt+1),&nz);
        if (nz != 0) {
/*  The second stage jumps directly to one of the third
 *  stage iterations and returns here if successful.
 *  Deflate the polynomial, store the zero or zeros and
 *  return to the main algorithm.
 */
            j = degree - n;
            
            zeror[j] = szr;
            
            zeroi[j] = szi;
            if( info != NULL )
               info[ j + 1 ] = itercnt;
            n -= nz;
            for (i=0;i<=n;i++)
                p[i] = qp[i];
                
            if (nz != 1) {
                zeror[j+1] = lzr;
                zeroi[j+1] = lzi; 
                if( info != NULL )
                  info[ j + 2 ] = 0;
            }
            goto _40;
        }
/*  If the iteration is unsuccessful another quadratic
 *  is chosen after restoring k.
 */
        for (i=0;i<n;i++) {
            k[i] = temp[i];
        }
    } 
/*  Return with failure if no convergence after 20 shifts. */
_99:
 
    delete [] svk;
    delete [] qk;
    delete [] k;
    delete [] qp;
    delete [] p;
    delete [] pt;
    delete [] temp;

   info[ 0 ] = clock() - sec;
   info[ 0 ] *= 1000;
   info[ 0 ] /= CLOCKS_PER_SEC;
 
    return degree - n;
}
/*  Computes up to L2 fixed shift k-polynomials,
 *  testing for convergence in the linear or quadratic
 *  case. Initiates one of the variable shift
 *  iterations and returns with the number of zeros
 *  found.
 */
void fxshfr(int l2,int *nz)
{
    double svu,svv,ui,vi,s;
    double betas,betav,oss,ovv,ss,vv,ts,tv;
    double ots,otv,tvv,tss;
    int type, i,j,iflag,vpass,spass,vtry,stry;

    *nz = 0;
    betav = 0.25;
    betas = 0.25;
    oss = sr;
    ovv = v;
/*  Evaluate polynomial by synthetic division. */
    quadsd(n,&u,&v,p,qp,&a,&b);
    calcsc(&type);
    for (j=0;j<l2;j++) {
/*  Calculate next k polynomial and estimate v. */
	nextk(&type);
	calcsc(&type);
	newest(type,&ui,&vi);
	vv = vi;
/*  Estimate s. */
        ss = 0.0;
        if (k[n-1] != 0.0) ss = -p[n]/k[n-1];
	tv = 1.0;
	ts = 1.0;
	if (j == 0 || type == 3) goto _70;
/*  Compute relative measures of convergence of s and v sequences. */
        if (vv != 0.0) tv = fabs((vv-ovv)/vv);
        if (ss != 0.0) ts = fabs((ss-oss)/ss);
/*  If decreasing, multiply two most recent convergence measures. */
	tvv = 1.0;
	if (tv < otv) tvv = tv*otv;
	tss = 1.0;
	if (ts < ots) tss = ts*ots;
/*  Compare with convergence criteria. */
	vpass = (tvv < betav);
	spass = (tss < betas);
	if (!(spass || vpass)) goto _70;
/*  At least one sequence has passed the convergence test.
 *  Store variables before iterating.
 */
	svu = u;
	svv = v;
	for (i=0;i<n;i++) {
		svk[i] = k[i];
	}
	s = ss;
/*  Choose iteration according to the fastest converging
 *  sequence.
 */
	vtry = 0;
	stry = 0;
	if (spass && (!vpass) || tss < tvv) goto _40;
_20:        
	quadit(&ui,&vi,nz);
        if (*nz > 0) return;
/*  Quadratic iteration has failed. Flag that it has
 *  been tried and decrease the convergence criterion.
 */
	vtry = 1;
	betav *= 0.25;
/*  Try linear iteration if it has not been tried and
 *  the S sequence is converging.
 */
	if (stry || !spass) goto _50;
	for (i=0;i<n;i++) {
		k[i] = svk[i];
	}
_40:
	realit(s,nz,&iflag);
	if (*nz > 0) return;
/*  Linear iteration has failed. Flag that it has been
 *  tried and decrease the convergence criterion.
 */
	stry = 1;
	betas *=0.25;
	if (iflag == 0) goto _50;
/*  If linear iteration signals an almost double real
 *  zero attempt quadratic iteration.
 */
	ui = -(s+s);
	vi = s*s;
	goto _20;
/*  Restore variables. */
_50:
	u = svu;
	v = svv;
	for (i=0;i<n;i++) {
		k[i] = svk[i];
	}
/*  Try quadratic iteration if it has not been tried
 *  and the V sequence is convergin.
 */
	if (vpass && !vtry) goto _20;
/*  Recompute QP and scalar values to continue the
 *  second stage.
 */
        quadsd(n,&u,&v,p,qp,&a,&b);
	calcsc(&type);
_70:
	ovv = vv;
	oss = ss;
	otv = tv;
	ots = ts;
    }
}
/*  Variable-shift k-polynomial iteration for a
 *  quadratic factor converges only if the zeros are
 *  equimodular or nearly so.
 *  uu, vv - coefficients of starting quadratic.
 *  nz - number of zeros found.
 */
void quadit(double *uu,double *vv,int *nz)
{
    double ui,vi;
    double mp,omp,ee,relstp,t,zm;
    int type,i,j,tried;

    *nz = 0;
    tried = 0;
    u = *uu;
    v = *vv;
    j = 0;
/*  Main loop. */
_10:    
   itercnt++;
    quad(1.0,u,v,&szr,&szi,&lzr,&lzi);
/*  Return if roots of the quadratic are real and not
 *  close to multiple or nearly equal and of opposite
 *  sign.
 */
    if (fabs(fabs(szr)-fabs(lzr)) > 0.01 * fabs(lzr)) return;
/*  Evaluate polynomial by quadratic synthetic division. */
    quadsd(n,&u,&v,p,qp,&a,&b);
    mp = fabs(a-szr*b) + fabs(szi*b);
/*  Compute a rigorous bound on the rounding error in
 *  evaluating p.
 */
    zm = sqrt(fabs(v));
    ee = 2.0*fabs(qp[0]);
    t = -szr*b;
    for (i=1;i<n;i++) {
       ee = ee*zm + fabs(qp[i]);
    }
    ee = ee*zm + fabs(a+t);
    ee *= (5.0 *mre + 4.0*are);
    ee = ee - (5.0*mre+2.0*are)*(fabs(a+t)+fabs(b)*zm);   
    ee = ee + 2.0*are*fabs(t);
/*  Iteration has converged sufficiently if the
 *  polynomial value is less than 20 times this bound.
 */
    if (mp <= 20.0*ee) {
        *nz = 2;
        return;
    }
    j++;
/*  Stop iteration after 20 steps. */
    if (j > 20) return;
    if (j < 2) goto _50;
    if (relstp > 0.01 || mp < omp || tried) goto _50;
/*  A cluster appears to be stalling the convergence.
 *  Five fixed shift steps are taken with a u,v close
 *  to the cluster.
 */
    if (relstp < eta) relstp = eta;
    relstp = sqrt(relstp);
    u = u - u*relstp;
    v = v + v*relstp;
    quadsd(n,&u,&v,p,qp,&a,&b);
    for (i=0;i<5;i++) {
	calcsc(&type);
	nextk(&type);
    }
    tried = 1;
    j = 0;
_50:
    omp = mp;
/*  Calculate next k polynomial and new u and v. */
    calcsc(&type);
    nextk(&type);
    calcsc(&type);
    newest(type,&ui,&vi);
/*  If vi is zero the iteration is not converging. */
    if (vi == 0.0) return;
    relstp = fabs((vi-v)/vi);
    u = ui;
    v = vi;
    goto _10;
}
/*  Variable-shift H polynomial iteration for a real zero.
 *  sss - starting iterate
 *  nz  - number of zeros found
 *  iflag - flag to indicate a pair of zeros near real axis.
 */
void realit(double sss, int *nz, int *iflag)
{
    double pv,kv,t,s;
    double ms,mp,omp,ee;
    int i,j;

    *nz = 0;
    s = sss;
    *iflag = 0;
    j = 0;
/*  Main loop */
    while (1) {
        itercnt++;
        pv = p[0];
/*  Evaluate p at s. */
        qp[0] = pv;
        for (i=1;i<=n;i++) {
            pv = pv*s + p[i];
            qp[i] = pv;
        }
        mp = fabs(pv);
/*  Compute a rigorous bound on the error in evaluating p. */
        ms = fabs(s);
        ee = (mre/(are+mre))*fabs(qp[0]);
        for (i=1;i<=n;i++) {
            ee = ee*ms + fabs(qp[i]);
        }
/*  Iteration has converged sufficiently if the polynomial
 *  value is less than 20 times this bound.
 */
        if (mp <= 20.0*((are+mre)*ee-mre*mp)) {
            *nz = 1;
            szr = s;
            szi = 0.0;   return ;                     // HVE return added
        }
        j++;
/*  Stop iteration after 10 steps. */
        if (j > 10) return;
        if (j < 2) goto _50;
        if (fabs(t) > 0.001*fabs(s-t) || mp < omp) goto _50;
/*  A cluster of zeros near the real axis has been
 *  encountered. Return with iflag set to initiate a
 *  quadratic iteration.
 */
        *iflag = 1;  sss =s;                          // HVE sss=s added
        return;
/*  Return if the polynomial value has increased significantly. */
_50:
        omp = mp;
/*  Compute t, the next polynomial, and the new iterate. */
        kv = k[0];
        qk[0] = kv;
        for (i=1;i<n;i++) {
            kv = kv*s + k[i];
            qk[i] = kv;
        }
        if (fabs(kv) <= fabs(k[n-1])*10.0*eta) {         // HVE n -> n-1
/*  Use unscaled form. */
            k[0] = 0.0;
            for (i=1;i<n;i++) {
                k[i] = qk[i-1];
            }
        }
        else {
/*  Use the scaled form of the recurrence if the value
 *  of k at s is nonzero.
 */
            t = -pv/kv;
            k[0] = qp[0];
            for (i=1;i<n;i++) {
                k[i] = t*qk[i-1] + qp[i];
            }
        }
        kv = k[0];
        for (i=1;i<n;i++) {
            kv = kv*s + k[i];
        }
        t = 0.0;
        if (fabs(kv) > fabs(k[n-1]*10.0*eta)) t = -pv/kv;
        s += t;
    }
}

/*  This routine calculates scalar quantities used to
 *  compute the next k polynomial and new estimates of
 *  the quadratic coefficients.
 *  type - integer variable set here indicating how the
 *  calculations are normalized to avoid overflow.
 */
void calcsc(int *type)
{
/*  Synthetic division of k by the quadratic 1,u,v */    
    quadsd(n-1,&u,&v,k,qk,&c,&d);
    if (fabs(c) > fabs(k[n-1]*100.0*eta)) goto _10;
    if (fabs(d) > fabs(k[n-2]*100.0*eta)) goto _10;
    *type = 3;
/*  Type=3 indicates the quadratic is almost a factor of k. */
    return;
_10:
    if (fabs(d) < fabs(c)) {
        *type = 1;
/*  Type=1 indicates that all formulas are divided by c. */   
        e = a/c;
        f = d/c;
        g = u*e;
        h = v*b;
        a3 = a*e + (h/c+g)*b;
        a1 = b - a*(d/c);
        a7 = a + g*d + h*f;
        return;
    }
    *type = 2;
/*  Type=2 indicates that all formulas are divided by d. */
	e = a/d;
	f = c/d;
	g = u*b;
	h = v*b;
	a3 = (a+g)*e + h*(b/d);
	a1 = b*f-a;
	a7 = (f+u)*a + h;
}
/*  Computes the next k polynomials using scalars 
 *  computed in calcsc.
 */
void nextk(int *type)
{
    double temp;
	int i;

    if (*type == 3) {
/*  Use unscaled form of the recurrence if type is 3. */
        k[0] = 0.0;
        k[1] = 0.0;
        for (i=2;i<n;i++) {
            k[i] = qk[i-2];
        }
        return;
    }
    temp = a;
    if (*type == 1) temp = b;
    if (fabs(a1) <= fabs(temp)*eta*10.0) {
/*  If a1 is nearly zero then use a special form of the
 *  recurrence.
 */
        k[0] = 0.0;
        k[1] = -a7*qp[0];
        for(i=2;i<n;i++) {
            k[i] = a3*qk[i-2] - a7*qp[i-1];
        }   return;           // HVE return added
    }
/*  Use scaled form of the recurrence. */
    a7 /= a1;
    a3 /= a1;
    k[0] = qp[0];
    k[1] = qp[1] - a7*qp[0];
    for (i=2;i<n;i++) {
	k[i] = a3*qk[i-2] - a7*qp[i-1] + qp[i];
    }
}
/*  Compute new estimates of the quadratic coefficients
 *  using the scalars computed in calcsc.
 */
void newest(int type,double *uu,double *vv)
{
    double a4,a5,b1,b2,c1,c2,c3,c4,temp;

/* Use formulas appropriate to setting of type. */
    if (type == 3) {
/*  If type=3 the quadratic is zeroed. */
        *uu = 0.0;
        *vv = 0.0;
        return;
    }
    if (type == 2) {
        a4 = (a+g)*f + h;
        a5 = (f+u)*c + v*d;
    }
    else {
        a4 = a + u*b +h*f;
        a5 = c + (u+v*f)*d;
    }
/*  Evaluate new quadratic coefficients. */
    b1 = -k[n-1]/p[n];
    b2 = -(k[n-2]+b1*p[n-1])/p[n];
    c1 = v*b2*a1;
    c2 = b1*a7;
    c3 = b1*b1*a3;
    c4 = c1 - c2 - c3;
    temp = a5 + b1*a4 - c4;
    if (temp == 0.0) {
        *uu = 0.0;
        *vv = 0.0;
        return;
    }
    *uu = u - (u*(c3+c2)+v*(b1*a1+b2*a7))/temp;
    *vv = v*(1.0+c4/temp);
    return;
}

/*  Divides p by the quadratic 1,u,v placing the quotient
 *  in q and the remainder in a,b.
 */
void quadsd(int nn,double *u,double *v,double *p,double *q,
    double *a,double *b)
{
    double c;
    int i;

    *b = p[0];
    q[0] = *b;
    *a = p[1] - (*b)*(*u);
    q[1] = *a;
    for (i=2;i<=nn;i++) {
        c = p[i] - (*a)*(*u) - (*b)*(*v);
	q[i] = c;
	*b = *a;
	*a = c;
    }	
}
/*  Calculate the zeros of the quadratic a*z^2 + b1*z + c.
 *  The quadratic formula, modified to avoid overflow, is used 
 *  to find the larger zero if the zeros are real and both
 *  are complex. The smaller real zero is found directly from 
 *  the product of the zeros c/a.
 */
void quad(double a,double b1,double c,double *sr,double *si,
        double *lr,double *li)
{
    double b,d,e;

    if (a == 0.0) {         /* less than two roots */
        if (b1 != 0.0) *sr = -c/b1;
        else  *sr = 0.0;
        *lr = 0.0;
        *si = 0.0;
        *li = 0.0;
        return;
    }
    if (c == 0.0) {         /* one real root, one zero root */
        *sr = 0.0;
	*lr = -b1/a;
        *si = 0.0;
        *li = 0.0;
        return;
    }
/* Compute discriminant avoiding overflow. */
    b = b1/2.0;
    if (fabs(b) < fabs(c)) { 
        if (c < 0.0) e = -a;
	else e = a;
        e = b*(b/fabs(c)) - e;
        d = sqrt(fabs(e))*sqrt(fabs(c));
    }
    else {
	e = 1.0 - (a/b)*(c/b);
        d = sqrt(fabs(e))*fabs(b);
    }
    if (e < 0.0) {      /* complex conjugate zeros */
	*sr = -b/a;
	*lr = *sr;
        *si = fabs(d/a);
	*li = -(*si);
    }
    else {
        if (b >= 0.0) d = -d;  /* real zeros. */
		*lr = (-b+d)/a;
                *sr = 0.0;
                if (*lr != 0.0) *sr = (c/ *lr)/a;
                *si = 0.0;
                *li = 0.0;
	}
}

