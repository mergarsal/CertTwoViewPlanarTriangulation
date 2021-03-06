# Certifiable two-view planar triangulation

**Manuscript submitted**

This repo contains the code for the solvers 
and optimality certificates 
for the two-view planar triangulation problem 
explained in [this paper](NONE) [1]. 






**Authors:** 
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718), 
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536)


**License:** [GPLv3](NONE)


If you use this code for your research, please cite:

```
NONE
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






