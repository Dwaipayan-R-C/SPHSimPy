#include "../headers/atom_structure.h"
#include "../headers/smoothing_kernel.h"
#include "../headers/xyz.h"

#include <gtest/gtest.h>

TEST(SmoothingKernel, AreaIntegrationTest)
{
    auto [positions]{read_xyz("../../xyz_output/kernel_test/kernel_test.xyz")};
    double smoothing_length{5.2};
    std::string particle_type[25];   
    
    Atoms atoms(positions);
    Smoothing_kernel smoothing_kernel(smoothing_length);
    double kernel_test_weights{0};
    double area{1 / pow(2.6, 2)};
    const double pi = 3.14159265358979323846;
    // Eigen::ArrayXd test_neighbor_kernel{atoms.positions.cols()};
    double kernel_value;
    std::ofstream outdata("../../xyz_output/kernel_test/kernel_data.out");
    for (int i = 0; i < 25; ++i)
    {
        particle_type[i]='X';
        kernel_test_weights += smoothing_kernel.calculate_spline_kernel(atoms.positions.col(12), atoms.positions.col(i));
    }
    
    for (int j{0}; j < 50; j++)
    {
        double kernel_data{0};
        for (int i = 0; i < 25; ++i)
        {
            kernel_value = smoothing_kernel.calculate_spline_kernel(atoms.positions.col(12), atoms.positions.col(i));
            kernel_data+=kernel_value/area;            
            atoms.velocities(2, i) = kernel_value;
        }
        outdata<<"[ "<<j+1<<" , "<<kernel_data<<" ],"<<std::endl;
        if(j==0)
        {
            write_xyz("../../xyz_output/kernel_test/kernel_test_out_"+to_string(j)+".xyz", atoms, particle_type);
        }
        write_xyz("../../xyz_output/kernel_test/kernel_test_out_"+to_string(j+1)+".xyz", atoms, particle_type);
        if(j<5)
        {
            atoms.positions(2, 12)-=.52;
            atoms.positions(1, 12)-=.52;
        }
        else if (j>=5 && j<25)
        {
            atoms.positions(2, 12)+=.52;           
        }
        else if (j>=25 && j<35)
        {
            atoms.positions(1, 12)+=.52; 
            atoms.positions(2, 12)-=.52;          
        }
        else{
            atoms.positions(2, 12)-=.52;
            atoms.positions(1, 12)-=.52;
        }
        
        
    }
    
    //    std::cout<<kernel_test_weights/area<<std::endl;
    EXPECT_NEAR(kernel_test_weights, area, 0.001);
}

TEST(SmoothingKernel, FirstDerivativeSummation)
{
    auto [positions]{read_xyz("../../xyz_output/kernel_test/kernel_test.xyz")};
    double smoothing_length{5.20000001};
    Atoms atoms(positions);
    Smoothing_kernel smoothing_kernel(smoothing_length);
    Eigen::Array2d first_derivative{0, 0};
    for (int i = 0; i < 25; ++i)
    {
        Eigen::Vector3d distance_vector{atoms.positions.col(12) - atoms.positions.col(i)};
        double r_norm{distance_vector.norm()};
        first_derivative += smoothing_kernel.calculate_kernel_derivative(atoms.positions.col(12), atoms.positions.col(i), r_norm);
        // std::cout << first_derivative(0) << std::endl;
        // std::cout << first_derivative(1) << std::endl;
    }

    EXPECT_NEAR(first_derivative(0), 0, 0.00005);
    EXPECT_NEAR(first_derivative(1), 0, 0.00005);
}

TEST(SmoothingKernel, OuterProductDerivative)
{
    auto [positions]{read_xyz("../../xyz_output/kernel_test/Kernel_test_derivative/kernel_test_derivative.xyz")};
    double smoothing_length{5.2};
    Atoms atoms(positions);
    Smoothing_kernel smoothing_kernel(smoothing_length);
    Eigen::Array2d first_derivative;
    Eigen::MatrixXd first_derivative_1(2,2);
    Eigen::MatrixXd first_derivative_outprod{2,2};
    first_derivative_outprod.setZero();
    Eigen::MatrixXd expected_mat{{-1.013/pow(smoothing_length/2,2),0},{0,-1.013/pow(smoothing_length/2,2)}};
    for (int i = 0; i < 9; ++i)
    {        
        Eigen::Vector3d distance_vector{atoms.positions.col(4) - atoms.positions.col(i)};
        Eigen::MatrixXd distance_2d_vector(2,2);
        distance_2d_vector(0)=distance_vector(1);
        distance_2d_vector(1)=distance_vector(2);
        double r_norm{distance_vector.norm()};
        first_derivative = smoothing_kernel.calculate_kernel_derivative(atoms.positions.col(4), atoms.positions.col(i), r_norm);
        first_derivative_1(0)=first_derivative(0);
        first_derivative_1(1)=first_derivative(1);
        first_derivative_outprod+=distance_2d_vector*first_derivative_1.transpose();    
        
    }

    std::cout <<"Analytical Expected: \n" <<expected_mat << std::endl;
    std::cout <<"Outer Product: \n" <<first_derivative_outprod << std::endl;

}