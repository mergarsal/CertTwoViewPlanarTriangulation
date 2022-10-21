# Certifiable two-view planar triangulation

**Manuscript submitted**

This repo contains the code for the solvers 
and optimality certificates 
for the two-view planar triangulation problem 
explained in [this paper](https://doi.org/10.1016/j.cviu.2022.103570). 






**Authors:** 
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)


**License:** [GPLv3](https://github.com/mergarsal/CertTwoViewPlanarTriangulation/blob/main/LICENSE)


If you use this code for your research, please cite:

```
@article{GARCIASALGUERO2022103570,
title = {Certifiable algorithms for the two-view planar triangulation problem},
journal = {Computer Vision and Image Understanding},
pages = {103570},
year = {2022},
issn = {1077-3142},
doi = {https://doi.org/10.1016/j.cviu.2022.103570},
url = {https://www.sciencedirect.com/science/article/pii/S1077314222001485},
author = {Mercedes Garcia-Salguero and Javier Gonzalez-Jimenez},
keywords = {Two-view planar triangulation, Primal–dual solver, Optimality certificate, Sufficient optimality condition},
}
```



## Build
```
git clone --recursive https://github.com/mergarsal/CertTwoViewPlanarTriangulation.git
cd CertTwoViewPlanarTriangulation

mkdir build & cd build 

cmake .. 

make -jX

```

The compiled examples should be inside the `bin` directory. 
Run for example: 
```
        ./bin/example_base
```
 
## Install 
In `build` folder: 
```
        sudo make install
```

We also provide the uninstall script: 
```
        sudo make uninstall
```



## How to use the library in your project

See the example in the folder `example_install` 
for the basic elements. 
       
1. In your CMakeLists.txt, add the dependences:
```
        find_package(OptTwoPlanarTriangulation   REQUIRED)
```

2. For your executable, add the library in 
```
        target_link_libraries(XXX OptTwoPlanarTriangulation)
```


## Test
We provide a file with the test parameters 
in *basic_params.txt*. 
The test reads the parameters in this file 
and perform the experiment. 






