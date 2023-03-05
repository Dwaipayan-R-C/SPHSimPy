#include <iostream>
#include <numeric>

#include "../headers/force_utils.h"


const void compute_density(NeighborList neighbor_list, Atoms &particles, double smoothing_length, int fluid_last_index)
{
    particles.densities.setZero();
    Smoothing_kernel smoothing_kernel(smoothing_length);
    for (auto[i, j]: neighbor_list)
    {
        particles.densities(i)+=particles.masses(j)*smoothing_kernel.calculate_spline_kernel(particles.positions.col(i),particles.positions.col(j));
    }
    particles.densities(Eigen::seq(fluid_last_index+1,Eigen::last)) = 1;
    
}

const void compute_pressure(NeighborList neighbor_list, Atoms &particles, double smoothing_length,int fluid_last_index,double rho_0, double K)
{
    particles.pressure=K*((particles.densities/rho_0)-1);
    for(int j{0};j<=fluid_last_index;j++)
    {
        if(particles.densities(j)<1)
        {
            particles.densities(j)=1;
        }
        if(particles.pressure(j)<0)
        {
            particles.pressure(j)=0;
        }
    }
}

const void compute_pressure_force(NeighborList neighbor_list, Atoms &particles, double smoothing_length, int fluid_last_index)
{
    Smoothing_kernel smoothing_kernel(smoothing_length);
    particles.pressure_forces.setZero();

    for (auto[i, j]: neighbor_list)
    {
        if(i<=fluid_last_index)
        {
            if(j<=fluid_last_index) {

                Eigen::Vector3d distance_vector = particles.positions.col(i) - particles.positions.col(j);
                double distance_{distance_vector.norm()};
                Eigen::Array2d W_del = smoothing_kernel.calculate_kernel_derivative(particles.positions.col(i),
                                                                                    particles.positions.col(j),
                                                                                    distance_);
                particles.pressure_forces(1, i) += -particles.masses(i) *
                                                   (particles.pressure(i) / pow(particles.densities(i), 2) +
                                                    particles.pressure(j) / pow(particles.densities(j), 2)) * W_del(0);
                particles.pressure_forces(2, i) += -particles.masses(i) *
                                                   (particles.pressure(i) / pow(particles.densities(i), 2) +
                                                    particles.pressure(j) / pow(particles.densities(j), 2)) * W_del(1);
            }
        }

    }
    for (auto[i, j]: neighbor_list)
    {
        if(i<=fluid_last_index)
        {
            if(j>fluid_last_index)
            {
                // particles.densities(j) = particles.densities(i);
                particles.pressure(j) = particles.pressure(i);
                Eigen::Vector3d distance_vector = particles.positions.col(i) - particles.positions.col(j);
                double distance_{distance_vector.norm()};
                Eigen::Array2d W_del=smoothing_kernel.calculate_kernel_derivative(particles.positions.col(i),particles.positions.col(j), distance_);
                particles.pressure_forces(1,i)+= -particles.masses(i)*(particles.pressure(i)/pow(particles.densities(i),2)+particles.pressure(j)/pow(particles.densities(j),2))*W_del(0);
                particles.pressure_forces(2,i)+= -particles.masses(i)*(particles.pressure(i)/pow(particles.densities(i),2)+particles.pressure(j)/pow(particles.densities(j),2))*W_del(1);

            }

        }

    }

}

const void compute_non_pressure_force(NeighborList neighbor_list, Atoms &particles, double smoothing_length,int fluid_last_index,double kin_viscosity_coeff)
{

    Smoothing_kernel smoothing_kernel(smoothing_length);
    particles.viscosity_forces.setZero();
    particles.velocities(Eigen::seq(0,2),Eigen::seq(fluid_last_index+1,Eigen::last))=0;
    for (auto[i, j]: neighbor_list)
    {
        if(i<=fluid_last_index)
        {
            Eigen::Vector3d distance_vector = particles.positions.col(i) - particles.positions.col(j);
            double distance_{distance_vector.norm()};
            Eigen::Array2d W_del=smoothing_kernel.calculate_kernel_derivative(particles.positions.col(i),particles.positions.col(j), distance_);
            Eigen::Array3d W_merge;
            W_merge(0)=0;
            W_merge(1)=W_del(0);
            W_merge(2)=W_del(1);

            if(j>fluid_last_index)
            {
                particles.densities(j)=particles.densities(i);
            }
            if(particles.positions(1,i)!=particles.positions(1,j) && particles.positions(2,i)!=particles.positions(2,j))
            {
                particles.viscosity_forces.col(i)+=W_merge*(particles.masses(j)/particles.densities(j)) * ((particles.velocities.col(i)-particles.velocities.col(j))*
                                                                                                           (particles.positions.col(i)-particles.positions.col(j)))/((particles.positions.col(i)-particles.positions.col(j))*(particles.positions.col(i)-particles.positions.col(j))+
                                                                                                                                                                     0.01*pow(2.6,2));
            }
//            std::cout<<particles.positions(1,j)<<std::endl;
//            std::cout<<particles.positions(2,j)<<std::endl;
//            std::cout<<particles.positions(1,i)<<std::endl;
//            std::cout<<particles.positions(2,i)<<std::endl;
//            std::cout<<particles.viscosity_forces(2,0)<<std::endl;
//            std::cout<<particles.densities(j)<<std::endl;
//            std::cout<<W_merge(1)<<std::endl;
//            std::cout<<W_merge(2)<<std::endl;
        }
    }
    particles.viscosity_forces(Eigen::seq(0,2),Eigen::seq(0,fluid_last_index)) *=2*kin_viscosity_coeff;

}

const void update_velocity_and_position(Atoms &particles, double del_t, Eigen::Array3Xd forces, int fluid_last_index)
{

    forces.row(0) = 0;
    particles.velocities(Eigen::seq(0,2),Eigen::seq(0,fluid_last_index)) += del_t*forces(Eigen::seq(0,2),Eigen::seq(0,fluid_last_index));
    particles.positions(Eigen::seq(0,2),Eigen::seq(0,fluid_last_index)) += del_t*particles.velocities(Eigen::seq(0,2),Eigen::seq(0,fluid_last_index));
}
