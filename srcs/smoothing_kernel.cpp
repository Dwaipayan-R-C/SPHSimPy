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
#include "../headers/smoothing_kernel.h"
#include <iostream>

Smoothing_kernel::Smoothing_kernel(double kernel_support) : kernel_support_{kernel_support}{}

const double Smoothing_kernel::calculate_spline_kernel(Eigen::Array3d particle_i, Eigen::Array3d particle_j)
{
    Eigen::Vector3d distance_vector = particle_i - particle_j;
    double distance_{distance_vector.norm()};    
    const double pi = 3.14159265358979323846;
    const double constant_coeff = 5 / (14 * pi * pow(kernel_support_/2, 2));
    
    // Kernel q calculation : q=||r||/h
    double kernel_weight{0};
    double q_{distance_ / (kernel_support_/2)};

    // Kernel conditions with smoothing length 
    if (0<=q_ && q_<1)
    {
        kernel_weight = constant_coeff*(pow((2-q_),3)-4*pow((1-q_),3));
    }
    else if (q_ >= 1 && q_ < 2)
    {
        kernel_weight= constant_coeff*pow((2-q_),3);
    }
    else{
        kernel_weight=0;
    }

    return kernel_weight;
}


const Eigen::Array2d Smoothing_kernel::calculate_kernel_derivative(Eigen::Array3d particle_i, Eigen::Array3d particle_j, double distance_norm)
{

    const double pi = 3.14159265358979323846;
    const double alpha = 5 / (14 * pi * pow(kernel_support_/2, 2));
    Eigen::Array2d displacement_vector{particle_i(1)-particle_j(1),particle_i(2)-particle_j(2)};
    Eigen::Array2d alpha_w = alpha*displacement_vector/(distance_norm*(kernel_support_/2));
    double q_{distance_norm / (kernel_support_/2)};
    Eigen::Array2d kernel_derivative{0,0};

    if(particle_i(1)==particle_j(1) && particle_i(2)==particle_j(2))
    {
        kernel_derivative(0)=0;
        kernel_derivative(1)=0;
    }
    else{
        if (0<=q_ && q_<1)
        {
            kernel_derivative = alpha_w*(-3*pow((2-q_),2)+12*pow((1-q_),2));
        }
        else if (q_ >= 1 && q_ < 2)
        {
            kernel_derivative= -alpha_w*3*pow((2-q_),2);
        }
    }
//    std::cout<<alpha_w(0)<<std::endl;
//    std::cout<<alpha_w(1)<<std::endl;
//    std::cout<<kernel_derivative(0)<<std::endl;
//    std::cout<<kernel_derivative(1)<<std::endl;
//    std::cout<<displacement_vector(0)<<std::endl;
//    std::cout<<displacement_vector(1)<<std::endl;
//    std::cout<<""<<std::endl;

    return kernel_derivative;
}



const double Smoothing_kernel::calculate_kernel_derivative2(Eigen::Array3d particle_i, Eigen::Array3d particle_j)
{
    Eigen::Vector3d distance_vector = particle_i - particle_j;
    double distance_{distance_vector.norm()};
    const double pi = 3.14159265358979323846;
    const double constant_coeff = 5 / (14 * pi * pow(kernel_support_, 2));

    // Kernel q calculation : q=||r||/h
    double kernel_weight{0};
    double q_{distance_ / (kernel_support_/2)};

    // Kernel conditions with smoothing length
    if (0<=q_ && q_<1)
    {
        kernel_weight = constant_coeff*(6*(2-q_)-24*(1-q_));
    }
    else if (q_ >= 1 && q_ < 2)
    {
        kernel_weight= constant_coeff*6*(2-q_);
    }
    else{
        kernel_weight=0;
    }
    return kernel_weight;
}