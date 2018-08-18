# cetsp-bezier-gsoa

This repository provides sources of the unsupervised learning based method for solving variants of the Traveling Salesman Problem (TSP) motivated by surveillance missions with aerial vehicles.
In particular, the problem is formulated as the TSP with the disk-shaped neighborhoods. The problem is further called the Close Enough TSP (CETSP) that is extended to consider the motion constraints of the multi-rotor aerial vehicles that are limited by the maximal velocity and acceleration.
Therefore, the requested trajectory to visit a given set of regions is considered as a sequence of Bézier curves. The proposed solver employs unsupervised learning that has been based initially on the Self-Organizing Map (SOM) but later consolidated and generalized to the Growing-Self-Organizing Array (GSOA) as a general approach for solving various routing problems.
The solver can solve 2D and also 3D instances where the shape of the neighborhood is represented by a sphere.
The approach has been introduced at ICRA 2018 and published in 

```
@article{faigl18ral,
   author    = {Jan Faigl and Petr Váňa},
   title     = {Surveillance Planning With B{\'{e}}zier Curves},
   journal   = {{IEEE} Robotics and Automation Letters},
   volume    = {3},
   number    = {2},
   pages     = {750--757},
   year      = {2018},
   url       = {https://doi.org/10.1109/LRA.2018.2789844},
   doi       = {10.1109/LRA.2018.2789844},
}
```

The approach has been further generalized for solving instances with a team of aerial vehicles; published in our article in the Journal of Field Robotics

```
@article{faigl18jfr,
   author    = {Jan Faigl and Petr Vana and Robert Penicka and Martin Saska},
   title     = {Unsupervised Learning based Flexible Framework for Surveillance Planning with Aerial Vehicles},
   journal   = {Journal of Field Robotics},
   year      = {2018},
   note      = {In press}
}
```

All the sources including the sub-merge crl can be compiled into the build directory using cmake by calling ./install-cmake_build.sh which creates the binary tcetsp-bezier-gsoa in the build directory.
Parameters of the solver can be provided from the command line or a configuration file, e.g., tcetsp-bezier-gsoa.cfg. The graphical output can be disabled by --gui none argument.

