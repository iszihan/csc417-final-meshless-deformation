
## CSC417 Final Project: Reimplementation of Meshless Deformation based on Shape Matching
## Instructions
To run our codes with different setups presented in the video, run the executable followed by an integer indicating the scenarios. 
- 0 for bunny collision
- 1 for cube collision 
- 2 for cubes with different clustering
- 3 for bunny animation while fixed to the floor

Example command:
  
    ./meshless 0

We allow for some interactivity:
* Press 1,2,or 3 to drop different objects into the scene.
* Press R to toggle to rigid body integration 
* Press L to toggle to linear integration. 
* Press Q to toggle to quadratic integration.

### Prerequisite installation

On all platforms, we will assume you have installed cmake and a modern c++
compiler on Mac OS X[¹](#¹macusers), Linux[²](#²linuxusers), or
Windows[³](#³windowsusers).

We also assume that you have cloned this repository using the `--recursive`
flag (if not then issue `git submodule update --init --recursive`). 

## Compilation for Debugging

This and all following assignments will follow a typical cmake/make build
routine. Starting in this directory, issue:

    mkdir build
    cd build
    cmake ..

If you are using Mac or Linux, then issue:

    make
