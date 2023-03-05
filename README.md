# Introduction
This project is developed for a course named - "Simulation in computer graphics" at University of Freiburg.
This skeleton uses [CMake](https://cmake.org/) and  automatically downloads [Eigen3](https://eigen.tuxfamily.org/) and
[Googletest](https://github.com/google/googletest).

To accomplish this course successfully, I cordially thank [Computer graphics department](https://cg.informatik.uni-freiburg.de/index.htm) to guide through the process. Also we use [Ovito](https://www.ovito.org/about/) to visualize our atoms and use  [Numpy](https://numpy.org/), [Matplotlib](https://matplotlib.org/) for plots.

# Visualization
Some plots are shown here 

<p align="center">
  <img src="figs/sph_3600particles.gif" alt="animated" width="200" />
</p>

Different timesteps of a 1 million particles simulation
<p align="center">
  <img src="figs/image_1.PNG" width="300"/>
</p>

## Run
To run the project, navigate to the source code and create this directory,
```
mkdir build                                                         # Create a build folder
cd build                                                            # inside that folder
cmake ..                                                            # configure the cmake filer
make                                                                # build the project
./myproject                                                         # run the make file
tests/myproject_tests
```  