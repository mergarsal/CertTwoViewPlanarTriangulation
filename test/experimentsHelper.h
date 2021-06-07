#pragma once 

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>  // for the file

#include <eigen3/Eigen/Dense>

using namespace std;


struct SceneOptions {
        double max_iter = 100;  // # Number of iterations 
        
        int n_size_img = 1;   // number of elements in FoV
        vector<double> size_imgs = {512};           // # FoV

        int n_parallax = 1; 
        vector<double> max_parallax = {2};  // # max parallax
        
        int n_dist_centers = 1; 
        vector<double> dist_centers = {4.0};  // # max parallax
        
        int n_focal = 1; 
        vector<double> focal_length = {1024}; 
        
        int n_noise = 1; 
        vector<double> noise = {0.5}; 
                
        bool use_custom_t = false; 
        
        Eigen::Matrix<double, 3, 1> custom_t; 
        
        double max_rotation = 0.5; 
        
        int method_trans = 1;
        
        
};


template <typename Scalar=double>
void readSetParams(ifstream & in_file, int & N, vector<Scalar> & arr);




bool readOptionsFromFile(string name_file, 
                         SceneOptions & options);




