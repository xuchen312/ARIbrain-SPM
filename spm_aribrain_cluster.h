#ifndef SPM_ARIBRAIN_CLUSTER_H
#define SPM_ARIBRAIN_CLUSTER_H

std::vector<double> findsimesfactor(bool, int);
std::vector<int> findhull(int, matlab::data::TypedArray<double>&);
std::vector<double> findalpha(matlab::data::TypedArray<double>&, int, std::vector<double>&, bool);
int findHalpha (std::vector<double>&, double, int);
int findConcentration(matlab::data::TypedArray<double>&, double, int, double, int);
int Find(int, std::vector<int>&);
void Union(int, int, std::vector<int>&, std::vector<int>&, std::vector<int>&);
int getCategory(double, double, double, int);
std::vector<int> findDiscoveries(matlab::data::TypedArray<int32_t>&, matlab::data::TypedArray<double>&, double, int, int, double, int, int);

#endif
