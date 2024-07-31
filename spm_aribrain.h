#ifndef SPM_ARIBRAIN_H
#define SPM_ARIBRAIN_H

int Find(int, std::vector<int>&);
void Union(int, int, std::vector<int>&, std::vector<int>&, std::vector<int>&);
int getCategory(double, double, double, int);
int findConcentration(std::vector<double>&, double, int, double, int);
std::vector<int> findDiscoveries(std::vector<int>&, std::vector<double>&, double, int, double, int, int);

#endif
