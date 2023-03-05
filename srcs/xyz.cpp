//
// Created by Dwaipayan on 15-05-2022.
//
#include <iomanip>
#include <iostream>
#include "../headers/xyz.h"

std::tuple<Positions_t> read_xyz(const std::string &filename) {
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;
        string name = "";
        // First line contains number of atoms
        std::getline(file, line);
        int nb_atoms;
        std::istringstream(line) >> nb_atoms;

        // Second line contains a comment - in extended XYZ this line contains
        // auxiliary information
        std::getline(file, line);

        // Data structures for names and positions
//        Names_t names(nb_atoms);
        Eigen::Array3Xd positions(3, nb_atoms);
        // Eigen::Array3Xd velocities(3, nb_atoms);
        positions.setZero();
        // velocities.setZero();

        // Now follows a line for each atom
        for (int i = 0; i < nb_atoms; ++i) {
            std::getline(file, line);
            std::istringstream(line) >>name>> positions(0, i) >> positions(1, i) >> positions(2, i);

        }

        // Close file, we're done
        file.close();

        return { positions};
    } else {
        throw std::runtime_error("Could not open file");
    }
}


std::tuple<Positions_t, Velocities_t> read_xyz_with_velocities(const std::string &filename) {
    std::ifstream file(filename);

    if (file.is_open()) {
        std::string line;
        string name = "";
        // First line contains number of atoms
        std::getline(file, line);
        int nb_atoms;
        std::istringstream(line) >> nb_atoms;

        // Second line contains a comment - in extended XYZ this line contains
        // auxiliary information
        std::getline(file, line);

        // Data structures for names and positions
//        Names_t names(nb_atoms);
        Eigen::Array3Xd positions(3, nb_atoms);
        Eigen::Array3Xd velocities(3, nb_atoms);
        positions.setZero();
        velocities.setZero();

        // Now follows a line for each atom
        for (int i = 0; i < nb_atoms; ++i) {
            std::getline(file, line);
            std::istringstream(line) >>name>> positions(0, i) >> positions(1, i) >> positions(2, i)
                                     >> velocities(0, i) >> velocities(1, i) >> velocities(2, i);

        }

        // Close file, we're done
        file.close();

        return { positions, velocities};
    } else {
        throw std::runtime_error("Could not open file");
    }
}

void write_xyz_vector(std::ofstream &file, Eigen::Array3Xd &vector_file) {
    // Number of atoms
    file << vector_file.cols() << std::endl;

    // Comment line
    file << std::endl;

    // Element name, position
    for (int i = 0; i < vector_file.cols(); ++i) {
        file << std::setw(2) << "X" << " "
             << std::setw(10) << vector_file.col(i).transpose()
             << std::endl;
    }
}

void write_xyz(std::ofstream &file, Atoms &atoms, string particle_type[]) {
    // Number of atoms
    file << atoms.nb_atoms() << std::endl;

    // Comment line
    file << std::endl;

    // Element name, position
    for (int i = 0; i < atoms.nb_atoms(); ++i) {
        file << std::setw(2) << particle_type[i] << " "
             << std::setw(10) << atoms.positions.col(i).transpose()
             << std::setw(10) << atoms.velocities.col(i).transpose()
             << std::endl;
    }
}


void write_xyz(const std::string &filename, Atoms &atoms, string particle_type[]) {
    // Open file for writing
    std::ofstream file(filename);

    // Write XYZ
    write_xyz(file, atoms,particle_type);

    // Close file, we're done
    file.close();
}


void write_xyz_vector(const std::string &filename, Eigen::Array3Xd &vector_file) {
    // Open file for writing
    std::ofstream file(filename);

    // Write XYZ
    write_xyz_vector(file, vector_file);

    // Close file, we're done
    file.close();
}
