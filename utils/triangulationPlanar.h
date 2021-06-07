#pragma once

#include <Eigen/Core>

namespace TwoViewPlanar
{
        typedef Eigen::Matrix<double, 2, 1> Vector2;
        typedef Eigen::Matrix<double, 3, 1> Vector3;
        typedef Eigen::Matrix<double, 3, 3> Matrix3;




        // TODO: move to utils
        void computeHomography(Vector3 & n, double d_plane, 
                               Matrix3 & R1, Matrix3 & R2, 
                               Vector3 & t1, Vector3 & t2, 
                               double focal_length, double size_img, 
                               Matrix3 & H, Matrix3 & Hp); 

        
        // Poly solver
        int correctMatchesPoly (const Matrix3& H,
                                const Vector2& point1,
                                const Vector2& point2,
                                Vector2* corrected_point1,
                                Vector2* corrected_point2);

        // solver Kanatani
        int correctMatchesKanatani (const Matrix3& H,
                                const Vector2& point1,
                                const Vector2& point2,
                                const double f0, 
                                Vector2* corrected_point1,
                                Vector2* corrected_point2);



         // general triangulation solver (polynomial)
         int FindObsHS (const Eigen::Matrix<double, 3, 3>& ematrix,
                                    const Eigen::Matrix<double, 2, 1>& point1,
                                    const Eigen::Matrix<double, 2, 1>& point2,
                                    Eigen::Matrix<double, 2, 1>* corrected_point1,
                                    Eigen::Matrix<double, 2, 1>* corrected_point2);
                            
                            
}  // end of namespace
