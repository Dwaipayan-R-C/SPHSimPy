/*
* Copyright 2023 Dwaipayan Roy Chowdhury
*
* ### MIT license
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/
#include <iostream>
#include "atom_structure.h"
#include "neighbors.h"
#include "smoothing_kernel.h"

#ifndef SPHPROJECT_FORCE_UTILS_H
#define SPHPROJECT_FORCE_UTILS_H

     /*
     * Calculates the density function  within the neighbor list
     */
const void compute_density(NeighborList neighbor_list, Atoms &particles,double smoothing_length, int fluid_last_index);

    /*
    * Calculates the density function  within the neighbor list
    */
const void compute_pressure(NeighborList neighbor_list, Atoms &particles, double smoothing_length,int fluid_last_index, double rho_0, double K = 2000);

    /*
    * Calculates the density function  within the neighbor list
    */
const void compute_pressure_force(NeighborList neighbor_list, Atoms &particles, double smoothing_length, int fluid_last_index);


/*
* Calculates the density function  within the neighbor list
*/
const void compute_non_pressure_force(NeighborList neighbor_list, Atoms &particles, double smoothing_length, int fluid_last_index,double kin_viscosity_coeff = 3.5);

const void update_velocity_and_position(Atoms &particles, double del_t,Eigen::Array3Xd forces, int fluid_last_index);


#endif //SPHPROJECT_FORCE_UTILS_H


