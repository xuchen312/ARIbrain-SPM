#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>


//------------------------- (0) PREPARATIONS FOR 3D INPUTS -------------------------//

// Macros:
// 1) Convert xyz coordinates to index
#define xyz2index(x, y, z, DIMS) ( (z-1)*DIMS[0]*DIMS[1] + (y-1)*DIMS[0] + (x-1) )
// 2) Compute size of 3D image
#define ndims(DIMS) ( DIMS[0]*DIMS[1]*DIMS[2] )

// Convert voxel index to [x y z] coordinates
std::vector<int> index2xyz(int               index,
                           std::vector<int>& DIMS)
{
    std::vector<int> XYZ(3);
    XYZ[0] = index % DIMS[0] + 1;
    //XYZ[0] = ((index-(XYZ[0]-1))/DIMS[0]) % DIMS[1] + 1;
    XYZ[1] = (index/DIMS[0]) % DIMS[1] + 1;
    //XYZ[1] = (index-(XYZ[0]-1)-(XYZ[1]-1)*DIMS[0]) / (DIMS[0]*DIMS[1]) + 1;
    XYZ[2] = index/(DIMS[0]*DIMS[1]) + 1;
    
    return XYZ;
}

// Convert several voxel indices to an xyz-coordinate matrix
Rcpp::IntegerMatrix ids2xyz(std::vector<int>& IDS,
                            std::vector<int>& DIMS)
{
    Rcpp::IntegerMatrix XYZS(IDS.size(), 3);
    for (int i = 0; i < IDS.size(); i++)
    {
        std::vector<int> XYZ = index2xyz(IDS[i], DIMS);
        XYZS(i,0) = XYZ[0];
        XYZS(i,1) = XYZ[1];
        XYZS(i,2) = XYZ[2];
    }

    return XYZS;
}

// Check if a voxel is in the mask
bool xyz_check(int               x,
               int               y,
               int               z,
               int               index,
               std::vector<int>& DIMS,
               std::vector<int>& MASK)
{
    return (x >= 1 && x <= DIMS[0] && \
            y >= 1 && y <= DIMS[1] && \
            z >= 1 && z <= DIMS[2] && \
            MASK[index] != 0);
}

// Find valid neighbours of a voxel
std::vector<int> findNeighbours(std::vector<int>& MASK,   // 3D mask of original orders (1:m)
                                std::vector<int>& DIMS,   // image dimensions
                                int               index,  // voxel index of interest
                                int               conn)   // connectivity criterion
{
    // compute [x y z] coordinates based on voxel index
    std::vector<int> XYZ = index2xyz(index, DIMS);
    
    // xyz coordinate adjustment vectors
    int DX[26] = {1,-1,0, 0,0, 0,  1,-1, 1,-1,1,-1, 1,-1,0, 0, 0, 0,  1,-1, 1,-1, 1,-1, 1,-1};
    int DY[26] = {0, 0,1,-1,0, 0,  1, 1,-1,-1,0, 0, 0, 0,1,-1, 1,-1,  1, 1,-1,-1, 1, 1,-1,-1};
    int DZ[26] = {0, 0,0, 0,1,-1,  0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1,  1, 1, 1, 1,-1,-1,-1,-1};
    
    // find all valid neighbours of a voxel
    std::vector<int> IDS;
    IDS.reserve(conn);
    for (int i = 0; i < conn; i++)
    {
        int      id = xyz2index(XYZ[0]+DX[i], XYZ[1]+DY[i], XYZ[2]+DZ[i], DIMS);
        bool inmask = xyz_check(XYZ[0]+DX[i], XYZ[1]+DY[i], XYZ[2]+DZ[i], id, DIMS, MASK);
        if (inmask)
        {
            IDS.push_back(MASK[id]);
        }
    }
    
    return IDS;
}

// Find the adjacency list for all in-mask voxels
std::vector<std::vector<int>> findAdjList(std::vector<int>& MASK,    // 3D mask of original orders (1:m)
                                          std::vector<int>& INDEXP,  // voxel indices of unsorted p-values
                                          std::vector<int>& DIMS,    // image dimensions
                                          int               m,       // number of in-mask voxels
                                          int               conn)    // connectivity criterion
{
    std::vector<std::vector<int>> ADJ;
    ADJ.reserve(m*conn);
    for (int i = 0; i < m; i++)
    {
        std::vector<int> IDS = findNeighbours(MASK, DIMS, INDEXP[i], conn);
        ADJ.push_back(IDS);
    }
    
    return ADJ;
}
