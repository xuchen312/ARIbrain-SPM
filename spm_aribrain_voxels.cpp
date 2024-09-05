#include "MatlabDataArray.hpp"
#include <vector>

#include "spm_aribrain_tdp.h"

//------------------------- PREPARATIONS FOR 3D INPUTS -------------------------//

// Convert xyz coordinates to 3D voxel index (starts from 0)
int xyz2index(int                                x,
              int                                y,
              int                                z,
              matlab::data::TypedArray<int32_t>& DIMS)  // 3D image dimensions
{
    return ( (z-1)*DIMS[0]*DIMS[1] + (y-1)*DIMS[0] + (x-1) );
}

// Convert 3D voxel index to [x y z] coordinates
std::vector<int> index2xyz(int                                index,  // 3D voxel index (starts from 0)
                           matlab::data::TypedArray<int32_t>& DIMS)   // 3D image dimensions
{
    std::vector<int> XYZ(3);
    XYZ[0] = index % DIMS[0] + 1;
    //XYZ[0] = ((index-(XYZ[0]-1))/DIMS[0]) % DIMS[1] + 1;
    XYZ[1] = (index/DIMS[0]) % DIMS[1] + 1;
    //XYZ[1] = (index-(XYZ[0]-1)-(XYZ[1]-1)*DIMS[0]) / (DIMS[0]*DIMS[1]) + 1;
    XYZ[2] = index/(DIMS[0]*DIMS[1]) + 1;
    
    return XYZ;
}

// Check if a voxel is in the mask
bool xyz_check(int                                x,
               int                                y,
               int                                z,
               int                                index,  // 3D voxel index (starts from 0)
               matlab::data::TypedArray<int32_t>& DIMS,   // 3D image dimensions
               matlab::data::TypedArray<int32_t>& MASK)   // 3D mask (non-zero for in-mask voxels)
{
    return (x >= 1 && x <= DIMS[0] && \
            y >= 1 && y <= DIMS[1] && \
            z >= 1 && z <= DIMS[2] && \
            MASK[index] != 0);
}

// Find valid neighbours of a voxel
std::vector<int> findNeighbours(matlab::data::TypedArray<int32_t>& MASK,   // 3D mask of original orders (1:m)
                                matlab::data::TypedArray<int32_t>& DIMS,   // image dimensions
                                int                                index,  // 3D voxel index (starts from 0)
                                int                                conn)   // connectivity criterion
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
            IDS.push_back(MASK[id]);  // append original/unsorted orders (1:m)
        }
    }
    
    return IDS;
}

// Find the adjacency list for all in-mask voxels
std::vector<std::vector<int>> findAdjList(matlab::data::TypedArray<int32_t>& MASK,    // 3D mask of original orders (1:m)
                                          matlab::data::TypedArray<int32_t>& INDEXP,  // in-mask voxel indices in 3D space (starts from 1)
                                          matlab::data::TypedArray<int32_t>& DIMS,    // image dimensions
                                          int                                m,       // number of in-mask voxels
                                          int                                conn)    // connectivity criterion
{
    std::vector<std::vector<int>> ADJ;
    ADJ.reserve(m*conn);
    for (int i = 0; i < m; i++)
    {
        std::vector<int> IDS = findNeighbours(MASK, DIMS, INDEXP[i]-1, conn);
        ADJ.push_back(IDS);
    }
    
    return ADJ;
}
