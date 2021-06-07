#include "generatePointCloudPlanar.h"


#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>


namespace TwoViewPlanar
{


Vector3 addNoise(double noise, double focal, size_t size_img, Vector3& obs);

Matrix3 generateRandomRotation( double maxAngle ); 

PCPlanarRes generatePointCloudPlanar(   PCPlanarParams & options, 
                                        const GenerateTranslation& generateTranslation, 
                                        const GenerateRotation& generateRotation,
                                        const PerturbTranslation& perturbTranslation, 
                                        const PerturbRotation& perturbRotation)
{
        
       
        size_t N_points = (options.N_width + 1) * (options.N_height + 1);
        PCPlanarRes res = PCPlanarRes(N_points);
        
        
        /** 1. GENERATE POINT CLOUD IN PLANE (N, D)**/
        
        double step_x = (options.width_plane / options.N_width * 2); 
        double step_y = (options.height_plane / options.N_height * 2); 
        
       
        
        res.points_3D.setZero();
        
        int idx_p = 0; 
        
        for (int i = 0; i <= options.N_width; i++)
        { 
                double px = step_x * i - options.width_plane; 
                
                for (int j = 0; j <= options.N_height; j++)
                { 
                        double py = step_y * j - options.height_plane; 
                        
                        // save point 
                        res.points_3D.col(idx_p) << px, py, options.d_plane; 
                        idx_p++;
                }  // end for Y
                
         }  // end for Xt
                
         
        
                
                
                /** 2. GENERATE POSE **/
                
                // 2.1 Pose for cam 1 (generic) 
                Vector3 t1 = generateRandomTranslationDefault(options.max_parallax, Vector3::Zero()); 
                Matrix3 R1 = generateRandomRotation(0.5); 
                        
                       
                
                // 2.2. Generate relative pose for second camera
               
                 
                Vector3 translation = generateTranslation(options.max_parallax, options.dir_parallax); 
                
                
                // std::cout << "[PC] Generate random rotation\n";  
                Matrix3 rotation = generateRotation(options.max_angle, options.dir_rotation); 
                

                
                // b. Perturb pose 
                translation = perturbTranslation(options.noise_trans, translation); 
                  
                rotation = perturbRotation(options.noise_rot, rotation); 
                
              
                
                Matrix3 R2 = rotation * R1; 
                Vector3 t2 = translation + R2*R1.transpose()*t1;
                
                // Save pose
                res.R1 = R1; 
                res.R2 = R2; 
                res.t1 = t1; 
                res.t2 = t2; 
                res.translation = translation; 
                res.rotation = rotation; 
                
                /** 3. GENERATE CORRESPONDENCES **/
                res.obs1.setZero(); 
                res.obs2.setZero(); 
                
               
                for (int i=0; i < N_points; i++)
                {
                       
                        Vector3 p1_3d = res.points_3D.col(i); 
                        
                        Vector3 obs1 = R1 * p1_3d + t1;  
                        
                        Vector3 obs2 = R2 * p1_3d + t2;
                                        
                        
                        // add noise 
                        obs1 = addNoise(options.noise, options.focal_length, options.size_img, obs1); 
                        obs2 = addNoise(options.noise, options.focal_length, options.size_img, obs2); 
                        
                        // save data 
                        res.obs1.col(i) = obs1; 
                        res.obs2.col(i) = obs2;                        
                }

               
        return res;
}  // end of function 









Vector3 addNoise(double noise, double focal, size_t size_img, Vector3& obs)
{
        Vector3 noisy_obs = obs; 
        
        noisy_obs(0) /= noisy_obs(2);
        noisy_obs(1) /= noisy_obs(2);
        noisy_obs(2) = 1;
        
        Matrix3 K = Matrix3::Identity(); 
        K(0, 0) = focal; 
        K(1, 1) = focal; 
        K(0, 2) = size_img / 2; 
        K(1, 2) = size_img / 2; 
        
        noisy_obs = K * noisy_obs; 
         
        
        noisy_obs(0) += (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * noise;
        noisy_obs(1) += (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * noise;

        
        
        return noisy_obs;
}


Vector3 generateRandomVector(double std_vec)
{
        Vector3 v;
        v.setZero(); 
        v(0) = (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * std_vec; 
        v(1) = (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * std_vec; 
        v(2) = (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * std_vec; 
        v.normalize();
        
        return v; 
}

Matrix3 generateRotationRodrigues(double max_angle, const Vector3 & dir_r)
{
                double angle = max_angle * (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
                
                Matrix3 rot; 
                Vector3 dir_rot = dir_r.normalized();
                Matrix3 T; 
                // cross matrix
                T << 0, -dir_rot(2), dir_rot(1), 
                     dir_rot(2), 0, -dir_rot(0), 
                     -dir_rot(1), dir_rot(0), 0;
                
                rot = Matrix3::Identity() + sin(angle) * T + (1 - cos(angle)) * T * T;
                       
                return rot;
}


// generate random translation with the specified norm
Vector3 generateRandomTranslationDefault( double max_parallax, const Vector3 & dir_parallax)
{
        return (max_parallax * generateRandomVector(1.0));

}
Matrix3 generateRandomRotation( double maxAngle )
        {
          
          Vector3 rpy;
          rpy[0] = ((double) std::rand())/ ((double) RAND_MAX);
          rpy[1] = ((double) std::rand())/ ((double) RAND_MAX);
          rpy[2] = ((double) std::rand())/ ((double) RAND_MAX);

          rpy[0] = maxAngle*2.0*(rpy[0]-0.5);
          rpy[1] = maxAngle*2.0*(rpy[1]-0.5);
          rpy[2] = maxAngle*2.0*(rpy[2]-0.5);

          Matrix3 R1;
          R1(0,0) = 1.0;
          R1(0,1) = 0.0;
          R1(0,2) = 0.0;
          R1(1,0) = 0.0;
          R1(1,1) = cos(rpy[0]);
          R1(1,2) = -sin(rpy[0]);
          R1(2,0) = 0.0;
          R1(2,1) = -R1(1,2);
          R1(2,2) = R1(1,1);

          Matrix3 R2;
          R2(0,0) = cos(rpy[1]);
          R2(0,1) = 0.0;
          R2(0,2) = sin(rpy[1]);
          R2(1,0) = 0.0;
          R2(1,1) = 1.0;
          R2(1,2) = 0.0;
          R2(2,0) = -R2(0,2);
          R2(2,1) = 0.0;
          R2(2,2) = R2(0,0);

          Matrix3 R3;
          R3(0,0) = cos(rpy[2]);
          R3(0,1) = -sin(rpy[2]);
          R3(0,2) = 0.0;
          R3(1,0) =-R3(0,1);
          R3(1,1) = R3(0,0);
          R3(1,2) = 0.0;
          R3(2,0) = 0.0;
          R3(2,1) = 0.0;
          R3(2,2) = 1.0;

          Matrix3 rotation = R3 * R2 * R1;

          rotation.col(0) = rotation.col(0) / rotation.col(0).norm();
          rotation.col(2) = rotation.col(0).cross(rotation.col(1));
          rotation.col(2) = rotation.col(2) / rotation.col(2).norm();
          rotation.col(1) = rotation.col(2).cross(rotation.col(0));
          rotation.col(1) = rotation.col(1) / rotation.col(1).norm();
          return rotation;
        }



// generate random rotation with maxAngle
Matrix3 generateRandomRotationDefault( double max_angle, const Vector3 & dir_rot)
{
                
        return (generateRandomRotation( max_angle ));         

} 

// generate orbital rotation with maxAngle
Matrix3 generateOrbitalRotation( double d, const Vector3 & trans)
{
        // here: 
        // d: distance to center point cloud 
        // trans: translation vector
        // compute angle
        Vector3 d_c; 
        d_c << 0, 0, d; 
        
        Vector3 r = d_c - trans; 
        
      
        double theta = r(2) / r.squaredNorm(); 
        
        // Construct matrix
        Matrix3 R = Matrix3::Identity(); 
        R(0, 0) = cos(theta); 
        R(0, 2) = sin(theta); 
        R(2, 0) = -sin(theta); 
        R(2, 2) = cos(theta); 
        
        return (R);         

} 


// generate random perturbation for translation 
Vector3 perturbRandomTranslationDefault( double noise, const Vector3 & trans)
{
        double n_trans = trans.norm(); 
        
        Vector3 res = trans; 
        // perturbation 
        res(0) +=  (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * noise;
        res(1) +=  (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * noise;
        res(2) +=  (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0 * noise;
        
        return (n_trans * res.normalized());
}
// generate random perturbation for translation 
Matrix3 perturbRandomRotationDefault( double noise, const Matrix3 & rot)
{
     Matrix3 noise_rot = generateRandomRotation(noise); 
     
     return (noise_rot * rot);
}


Vector3 generateTranslationForward( double max_parallax, const Vector3 & dir_parallax)
{
        Vector3 t; 
        t << 0, 0, 1; 
        
        return (max_parallax * t);

}; 

Vector3 generateTranslationStereo( double max_parallax, const Vector3 & dir_parallax)
{
        Vector3 t; 
        t << 1, 0, 0; 
        
        return (max_parallax * t);

}; 


Vector3 generateTranslationSideways( double max_parallax, const Vector3 & dir_parallax)
{
        Vector3 t; 
        t(0) = (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0; 
        t(1) = (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0; 
        t(2) = 0; 
        t.normalize(); 
        
        return (max_parallax * t);

}; 

Vector3 generateTranslationOblique( double max_parallax, const Vector3 & dir_parallax)
{
        Vector3 t;
        t << 1, 1, 1; 
        t.normalize(); 
        
        return (max_parallax * t);

}; 

}  // end of namespace UtilsTwoView
