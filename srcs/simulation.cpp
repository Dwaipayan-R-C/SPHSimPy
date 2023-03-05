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

#include "../headers/simulation.h"
#include "../headers/neighbors.h"
#include "../headers/xyz.h"
#include "Eigen/Core"
#include <tuple>
#include <filesystem>
using std::ofstream;

// This function runs the main simulation
void sph_simulation(double particle_size, double del_t, int n_steps, int save_every, int fluid_last_index, double rho_0,double stiffness_K, double visc_coeff)
{
    // Input output files
    std::ofstream outdata("../xyz_output/data/"+to_string((int)(del_t*1000))+"_"+to_string((int)stiffness_K)+"_"+to_string((int)visc_coeff)+".out");
    if (!outdata) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
    }
    auto [positions]{read_xyz("../xyz_output/particles.xyz")};
    double rc = particle_size*2+0.000001;
    Atoms atoms(positions);    
    std::string particle_type[atoms.nb_atoms()];
    for (int j{0};j<=fluid_last_index;j++)
    {
        particle_type[j]="F";
    }
    for (int j{fluid_last_index+1};j<atoms.nb_atoms();j++)
    {
        particle_type[j]="B";
    }
    
    NeighborList neighbor_list(rc);
    Eigen::Array3d gravity{0,0,-9.81};

    std::cout<<"Pressure"<<" , "<<"Densities"<<" , "<<"Pressure forces"<<" , "<<"Viscous forces"<<" , "<<"Total forces"<<std::endl;
    int count{0};

    //Notions: particle dimension = 2 mm, density of water = 1000kgs/m3        
    atoms.masses = pow(particle_size,2)*rho_0;
    for(int i{0};i<n_steps;i++)
    {
        neighbor_list.update(atoms);
        compute_density(neighbor_list,atoms,rc,fluid_last_index);
        compute_pressure(neighbor_list,atoms,rc,fluid_last_index,rho_0, stiffness_K);
        compute_pressure_force(neighbor_list,atoms,rc, fluid_last_index);
        compute_non_pressure_force(neighbor_list,atoms,rc, fluid_last_index, visc_coeff);
        Eigen::Array3Xd total_force = atoms.pressure_forces + atoms.viscosity_forces;
        total_force(Eigen::seq(0,2),Eigen::seq(0,fluid_last_index)).row(2)-=9.81;
        std::cout<<i*del_t<<" , "<<atoms.pressure(0)<<" , "<<atoms.densities(0)<<" , "<<atoms.pressure_forces(2,0)<<" , "<<atoms.viscosity_forces(2,0)<<" , "<<total_force(2,0)<<std::endl;
        update_velocity_and_position(atoms,del_t,total_force, fluid_last_index);
        double total_density = atoms.densities(Eigen::seq(0,fluid_last_index)).sum()/atoms.densities(Eigen::seq(0,fluid_last_index)).rows();        
        double time = i*del_t;
        namespace fs = std::filesystem;
        fs::create_directories("../xyz_output/sph_output/"+to_string((int)(del_t*1000))+"_"+to_string((int)stiffness_K)+"_"+to_string((int)visc_coeff));
        outdata<<time<<" , "<<total_density<<std::endl;            
        if(i%save_every==0)
        {
            write_xyz("../xyz_output/sph_output/"+to_string((int)(del_t*1000))+"_"+to_string((int)stiffness_K)+"_"+to_string((int)visc_coeff)+"/sph_vector_"+to_string(count+1)+".xyz", atoms, particle_type);
            count+=1;
        }

    }

}