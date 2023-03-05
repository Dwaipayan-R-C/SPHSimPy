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
#include "headers/simulation.h"


int main() {
    double del_t_list[3]{0.008,0.01,0.02};    
    double stiffness_K_list[1]{30000};
    double visc_coeff_list[1]{18};
    int n_steps=3000;
    double particle_size = 2.6;
    int save_every=25;
    int fluid_last_index=399;
    for (int j=0;j<3;j++)
    {
        for (int l=0;l<1;l++)
        {
            for (int k=0;k<1;k++) 
            {                
                double del_t = del_t_list[j];   
                double stiffness_K=stiffness_K_list[l];
                double visc_coeff=visc_coeff_list[k];         
                double rho_0=1;                 // in kg/m3    
                sph_simulation(particle_size, del_t,  n_steps,  save_every,  fluid_last_index, rho_0,stiffness_K, visc_coeff);
            }
        }
    }
    return 0;
}
