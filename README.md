# DGGML: The Dynamical Graph Grammar Modeling Library
DGGML is a modeling library for developing simulations of spatially embedded dynamical graph grammars. More details of the algoirthm can be found in the arXiv version of the article, ["Advances in the Simulation and Modeling of Complex Systems using Dynamical Graph Grammars"](https://arxiv.org/abs/2407.10072). It's a follow up to [CajeteCMA](https://github.com/emedwede/CajeteCMA) - a Graph Grammar Simulator for the Cortical Microtubule Array (CMA), which has a serial implementation for a graph grammar simulator featured in the paper titled, ["Approximate Simulation of Cortical Microtubule Models using Dynamical Graph Grammars"](https://dx.doi.org/10.1088/1478-3975/acdbfb). The graph library used is [YAGL: Yet another Graph Library](https://github.com/emedwede/YAGL.git). Additional info can be found on my [thesis website](https://emedwede.github.io).

Dynamical Graph Grammars (DGGs) have a long history, and there are many papers discussing the theory. Lots of resources can be found [here](https://emj.ics.uci.edu), and a more recent open acess paper, ["Prospects for Declarative Mathematical Modeling of Complex Biological Systems"](https://doi.org/10.1007/s11538-019-00628-7) is a great place to start. Other realted computational biology software can be found [here](https://emj.ics.uci.edu/software/) and [Plenum](https://computableplant.ics.uci.edu/theses/guy/Plenum.html) is also worth checking out. It is the first implementation of DGGs, but requires Mathematica to use.  


# Requirements and Build Instructions
YAGL is a built in dependency, so please use `git clone https://github.com/emedwede/DGGML.git --recursive`

## General
The minimum CXX version required to compile is C++17 and to build CMake version 3.16 and above. The code has been tested on Linux and the latest version of MacOS using gcc12, but should be compilable on Windows as well. We also use the state of the art ODE solver [SUNDIALS](https://github.com/LLNL/sundials), and currently a build of version 6.7.0 is required, and will be downloaded an installed automatically if you do not have the required version. For the first install, this may take a few minutes to download. 

## MacOS
To get the code to run on the new ARM based Apple silicon(poorly named M1/M2/pro/max/etc), a few commands in homebrew will get you going. It'll probably work on older intel based chips too. If you haven't already, check out [homebrew](https://brew.sh/) and follow these instructions:

**Warining, brew may update your old packages when installing new ones, so be careful and read the docs if you're new**

1) brew install brew install gcc@12
2) brew install cmake 3.16 or greater
4) clone this repository to wherever you want using `git clone https://github.com/emedwede/DGGML.git --recursive`
6) cd into the repo, type `mkdir build && cd build` or wherever you choose to build.
7) from the build directory set the compilers and configure the build with `CC=gcc-12 CXX=g++12 cmake ..`
   1) Alternatively: cmake   `-D CMAKE_CXX_FLAGS="-O3” -D CMAKE_CXX_EXTENSIONS=Off  -D CMAKE_INSTALL_PREFIX=[Your/install/path] -D SUNDIALS_DIR=[Your/SUNDIALS/build/path [Your/Path]/DGGML`
8) Then to build `make -j4`
9) This will take a few minutes, and after complete you can navigate to `examples/CMA` and run `./mt_dgg_simulator settings.json` to see the example code in action. 

Please not that the clang compiler is currently untested, so use at your discretion.

## Linux
Nearly identical to the MacOS instructions. Instead build what you need from source or use your favorite package manager.

**Warining, brew may update your old packages when installing new ones, so be careful and read the docs if you're new**

1) install gcc12
2) install cmake 3.16 or greater
4) clone this repository to wherever you want using `git clone https://github.com/emedwede/DGGML.git --recursive`
6) cd into the repo, type `mkdir build && cd build` or wherever you choose to build.
7) from the build directory set the compilers and configure the build with `CC=gcc-12 CXX=g++12 cmake ..`
   1) Alternatively: cmake   `-D CMAKE_CXX_FLAGS="-O3” -D CMAKE_CXX_EXTENSIONS=Off  -D CMAKE_INSTALL_PREFIX=[Your/install/path] -D SUNDIALS_DIR=[Your/SUNDIALS/build/path [Your/Path]/DGGML`
8) Then to build `make -j4`
9) This will take a few minutes, and after complete you can navigate to `examples/CMA` and run `./mt_dgg_simulator settings.json` to see the example code in action.

# Usage
The library is meant to allow the user to build simulations of DGGs. Please see the CMA example for how to write grammar rules and build a simulator. Visualization, checkpointing, metric collection and most pieces are a choose your own adventure. You can use the defaults or write your own. I recommend using [Paraview](https://www.paraview.org) to visualize simulation results until you get the hang of it. The code includes a paraview VTK file writer by default, but this can also be overriden and you can also write your own.  
