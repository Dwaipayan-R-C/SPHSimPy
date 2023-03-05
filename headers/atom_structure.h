/*
* Copyright 2022 Dwaipayan Roy Chowdhury
*
* ### MIT license
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
#include "Eigen/Core"

#ifndef YAMD_ATOM_STRUCTURE_H
#define YAMD_ATOM_STRUCTURE_H

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Pressure_Force_t = Eigen::Array3Xd;
using Viscosity_Force_t = Eigen::Array3Xd;
using Masses_t = Eigen::ArrayXd;
using Pressure_t = Eigen::ArrayXd;
using Density_t = Eigen::ArrayXd;


using namespace std;

struct Atoms {
    Masses_t masses;
    Positions_t positions;
    Velocities_t velocities;
    Pressure_Force_t pressure_forces;
    Viscosity_Force_t viscosity_forces;
    Pressure_t pressure;
    Density_t densities;


    Atoms(Positions_t &p)
        : positions{p},
          velocities{3, p.cols()},
          pressure_forces{3, p.cols()},
          viscosity_forces{3, p.cols()},
          masses{p.cols()},
          pressure{p.cols()},
          densities{p.cols()} {
        velocities.setZero();
        pressure_forces.setZero();
        masses.setOnes();
        pressure.setZero();
        viscosity_forces.setZero();
        densities.setZero();
    }    
    
    Eigen::Index nb_atoms() const {
        return positions.cols();
    }
};
#endif // YAMD_ATOM_STRUCTURE_H
