#ifndef SPM_ARIBRAIN_TDP_H
#define SPM_ARIBRAIN_TDP_H

int xyz2index(int, int, int, matlab::data::TypedArray<int32_t>&);
std::vector<int> index2xyz(int, matlab::data::TypedArray<int32_t>&);
bool xyz_check(int, int, int, int, matlab::data::TypedArray<int32_t>&, matlab::data::TypedArray<int32_t>&);
std::vector<int> findNeighbours(matlab::data::TypedArray<int32_t>&, matlab::data::TypedArray<int32_t>&, int, int);
std::vector<std::vector<int>> findAdjList(matlab::data::TypedArray<int32_t>&, matlab::data::TypedArray<int32_t>&, matlab::data::TypedArray<int32_t>&, int, int);

#endif
