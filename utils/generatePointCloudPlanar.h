#pragma once 

#include <Eigen/Core>
#include <functional>


namespace TwoViewPlanar
{

typedef Eigen::Matrix<double, 2, 1> Vector2;
typedef Eigen::Matrix<double, 3, 1> Vector3;
typedef Eigen::Matrix<double, 3, 3> Matrix3;


/* Generate random translation */
using GenerateTranslation = std::function<Vector3(const double max_parallax, const Vector3& direction_parallax)>;
/* Generate random rotation */
using GenerateRotation = std::function<Matrix3(const double max_angle, const Vector3& rot_angle)>;
/* Generate random perturbation for translation */
using PerturbTranslation = std::function<Vector3(const double max_parallax, const Vector3& direction_parallax)>;
/* Generate random perturbation for rotation */
using PerturbRotation = std::function<Matrix3(const double max_angle, const Matrix3& rot)>;                
                               
// generate random translation with the specified norm
Vector3 generateRandomTranslationDefault( double max_parallax, const Vector3 & dir_parallax);
// generate random rotation with maxAngle
Matrix3 generateRandomRotationDefault( double maxAngle, const Vector3 & dir_rot); 


// generate random perturbation for translation 
Vector3 perturbRandomTranslationDefault( double noise, const Vector3 & trans);
// generate random perturbation for translation 
Matrix3 perturbRandomRotationDefault( double noise, const Matrix3 & rot); 

            

/* Options for the point cloud generation */
struct PCPlanarParams
{
                   EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                   
                   double focal_length = 800;           // in pixels
                   size_t size_img = 1024;
                   size_t N_height = 10;
                   size_t N_width = 10;
                   double noise = 0.0;                  // in pixels
                   double width_plane = 2;              // this is half the plane so: [-width, +width]
                   double height_plane = 2;    
                   
                   double d_plane = 3; 
                   
                   
                   // params for relpose
                   double max_parallax = 2.0;               // in meters
                   double max_angle = 0.5;           // in degrees
                   
                   double noise_trans = 0.0; 
                   double noise_rot = 0.0; 
                   
                   Vector3 dir_parallax; 
                   Vector3 dir_rotation;
                   
                   // constructor
                   PCPlanarParams(){}; 



};  // end of struct PCParams


/* Result of the point cloud generation */
struct PCPlanarRes
{
                 EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                   
                    // Poses
                    Vector3 t1; 
                    Vector3 t2;
                    Vector3 translation;
                    Matrix3 R1; 
                    Matrix3 R2; 
                    Matrix3 rotation;
                    
                    
                    // Observations
                    Eigen::MatrixXd  obs1;  // frame 0
                    Eigen::MatrixXd  obs2;  // frame 1
                    // 3D world points
                    Eigen::MatrixXd  points_3D;
                    
                    // homographies
                    Matrix3 H;   // normalized K = eye(3)
                    Matrix3 Hp;  // K
                    
                    // constructor
                    /// default
                    PCPlanarRes()
                    {
                    t1.setZero(); 
                    t2.setZero(); 
                    R1.setZero(); 
                    R2.setZero();
                    translation.setZero(); 
                    rotation.setZero(); 
                    points_3D = Eigen::MatrixXd(3, 10); 
                    obs1 = Eigen::MatrixXd(3, 10); 
                    obs2 = Eigen::MatrixXd(3, 10);  
                    
                    H.setZero(); 
                    Hp.setZero();                   
                    }; 
                     
                    PCPlanarRes(size_t n_points) {
                    
                    t1.setZero(); 
                    t2.setZero(); 
                    R1.setZero(); 
                    R2.setZero();
                    translation.setZero(); 
                    rotation.setZero();
                     
                    points_3D = Eigen::MatrixXd(3, n_points); 
                    obs1 = Eigen::MatrixXd(3, n_points); 
                    obs2 = Eigen::MatrixXd(3, n_points); 
                    
                    H.setZero(); 
                    Hp.setZero(); 
                    }; 
        
}; // end of PCRes  
  
  
PCPlanarRes generatePointCloudPlanar(PCPlanarParams & options, 
                         const GenerateTranslation& generateTranslation = generateRandomTranslationDefault,
                         const GenerateRotation& generateRotation = generateRandomRotationDefault,
                         const PerturbTranslation& perturbTranslation = perturbRandomTranslationDefault,
                         const PerturbRotation& perturbRotation = perturbRandomRotationDefault);
                         
 




/** Some special functions  **/
Vector3 generateTranslationForward( double max_parallax, const Vector3 & dir_parallax); 

Vector3 generateTranslationStereo( double max_parallax, const Vector3 & dir_parallax);

Vector3 generateTranslationSideways( double max_parallax, const Vector3 & dir_parallax); 

Vector3 generateTranslationOblique( double max_parallax, const Vector3 & dir_parallax); 

Vector3 generateTranslationOrbital( double max_parallax, const Vector3 & dir_parallax); 

Matrix3 generateOrbitalRotation( double max_parallax, const Vector3 & dir_parallax); 
                        
};  // end of namespace UtilsTwoView

