#include "../headers/atom_structure.h"
#include "../headers/neighbors.h"
#include "../headers/xyz.h"

#include <gtest/gtest.h>

TEST(NeighborsTest, NeighborListTest) {    
    // Positions_t positions(3, 1);
    
    auto [positions]{read_xyz("../../xyz_output/kernel_test/kernel_test.xyz")};
    double rc = 5.20000001;
    Atoms atoms(positions);
    NeighborList neighbor_list(rc);
    auto &[seed, neighbors]{neighbor_list.update(atoms)};

    // All atoms except 3 and 4 are neighbors of each other
    // EXPECT_EQ(neighbor_list.nb_neighbors(), 10);
    EXPECT_EQ(neighbor_list.nb_neighbors(12), 13);    // Center atom
    EXPECT_EQ(neighbor_list.nb_neighbors(12), 13);     // Last corner atom

    // EXPECT_TRUE((neighbors(Eigen::seq(seed(0), seed(1) - 1)) == Eigen::Array3i{3, 1, 2}).all());
    
}