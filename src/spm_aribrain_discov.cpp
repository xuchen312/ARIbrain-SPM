#include "MatlabDataArray.hpp"
#include <vector>
#include <cmath>

#include "spm_aribrain_cluster.h"

//------------------------- FUNCTIONS USED IN BOTH MEX FILES -------------------------//

// Calculate the denominator of the local test
std::vector<double> findsimesfactor(bool simes,  // assume simes yes or no
                                    int  m)      // number of p-values
{
  std::vector<double> simesfactors(m+1);  // denominator of simes test
  double multiplier = 0;                  // extra multiplier term needed for non-simes
  simesfactors[0] = 0;
  if (simes)
    for (int i=1; i<=m; i++)
      simesfactors[i] = i;
  else
    for (int i=1; i<=m; i++)
    {
      multiplier += 1.0/i;
      simesfactors[i] = i*multiplier;
    }

  return simesfactors;
}

// Implementation of Fortune 1989
std::vector<int> findhull(int                               m,      // length of p
                          matlab::data::TypedArray<double>& sortp)  // p-values (sorted!)
{
  // intialize output length
  int r;
  std::vector<int> hull(1);
  hull.push_back(1);
  bool notconvex;
  
  // find the hull
  for (int i=2; i<=m; i++)
  {
    if (i==m || (m-1)*(sortp[i-1]-sortp[0]) < (i-1)*(sortp[m-1]-sortp[0]))
    {
      do {
        r = hull.size()-1;
        if (r>1)
        {
          notconvex = ((i-hull[r-1])*(sortp[hull[r]-1]-sortp[hull[r-1]-1]) >= (hull[r]-hull[r-1])*(sortp[i-1]-sortp[hull[r-1]-1]));
        } else {
          if (r==1)
          {
            notconvex = (i*sortp[hull[1]-1] >= hull[1]*sortp[i-1]);
          } else {
            notconvex = false;
          }
        }
        if (notconvex) hull.resize(r);
      } while (notconvex);
      hull.push_back(i);
    }
  }
  
  return hull;
}

// Find the jumps of h(alpha)
std::vector<double> findalpha(matlab::data::TypedArray<double>& sortp,         // p-values (sorted!)
                              int                               m,             // length of p
                              std::vector<double>&              simesfactors,  // denominators of local test
                              bool                              simes)         // assume simes yes or no
{
  // initialize output
  std::vector<double> alpha(m+1);
  
  // initialize
  std::vector<int> hull = findhull(m,sortp);
  double Dk = 0;
  int     k = hull.size()-1;
  int     i = 1;
  
  // algorithm for alpha*
  while (i<=m)
  {
    if (k>1)
      Dk = sortp[hull[k-1]-1]*(hull[k]-m+i)-sortp[hull[k]-1]*(hull[k-1]-m+i);
    if (k>1 && Dk<0)
    {
      k--;
    } else {
      alpha[i-1] = simesfactors[i]*sortp[hull[k]-1]/(hull[k]-m+i);
      i++;
    }
  }
  
  // Bound alpha by 1
  if (!simes)
  {
    for (int i=m-1; i>=0; i--)
      if (alpha[i]>1)
        alpha[i] = 1;
  }
  
  // Find the cumulative maximum
  if (!simes)
  {
    for (int i=m-2; i>=0; i--) {
      if (alpha[i] < alpha[i+1])
        alpha[i] = alpha[i+1];
    }
  }
  
  // add (m+1)st element to alpha (Not needed as it's already added!)
  //alpha.push_back(0);
  
  return alpha;
}

// Calculate the value of h(alpha) for a given alpha
int findHalpha (std::vector<double>& jumpalpha,  // points where h jumps
                double               alpha,      // alpha where h is to be evaluated
                int                  m)          // size of the multiple testing problem
{
  // We use bisection
  int lower = 0;
  int upper = m+1;
  int mid   = 0;
  while (lower+1 < upper)
  {
    mid = (lower+upper+1)/2;
    if (jumpalpha[mid-1] > alpha)
    {
      lower = mid;
    }
    else
    {
      upper = mid;
    }
  }
  return lower;
}
 
// Calculates the size of the concentration set at a fixed alpha
int findConcentration(matlab::data::TypedArray<double>& sortp,   // p-values (sorted!)
                      double                            simesh,  // simesfactor at h(alpha)
                      int                               h,       // h(alpha)
                      double                            alpha,   // alpha itself
                      int                               m)       // size of the problem
  
{
  // from m-h we increase z until we fulfil the condition
  int z = m-h;
  if (z>0)  // h=m implies z=0
  {
    while ((z<m) & (simesh*sortp[z-1] > (z-m+h+1)*alpha))
    {
      z++;
    }
  }
  return z;
}

// Find function for disjoint set data structure
// (1) Old recursive version
// int Find(int               x,
//          std::vector<int>& parent)
// {
//   if (parent[x] != x)
//   {
//     parent[x] = Find(parent[x], parent);
//   }
//
//   return parent[x];
// }
// (2) iterative version (more stable)
int Find(int               x,
         std::vector<int>& parent)
{
  while (parent[x] != x)
  {
    parent[x] = parent[parent[x]];
    x         = parent[x];
  }
  
  return x;
}

// Union function for disjoint set data structure
// Extra: we keep track of the lowest entry of each disjoint set
// That way we can find the lower set to merge with
void Union(int               x,
           int               y,
           std::vector<int>& parent,
           std::vector<int>& lowest,
           std::vector<int>& rank)
{
  int xRoot = Find(x, parent);
  int yRoot = Find(y, parent);
  
  // if x and y are already in the same set (i.e., have the same root or representative)
  if (xRoot == yRoot) return; // Note: this never happens in our case
  
  // x and y are not in same set, so we merge
  if (rank[xRoot] < rank[yRoot])
  {
    parent[xRoot] = yRoot;
    lowest[yRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
  else if (rank[xRoot] > rank[yRoot])
  {
    parent[yRoot] = xRoot;
    lowest[xRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
  else
  {
    parent[yRoot] = xRoot;
    rank[xRoot]++;
    lowest[xRoot] = std::min(lowest[xRoot], lowest[yRoot]);
  }
}

// Calculate the category for each p-value
int getCategory(double p,       // p-value for which we need the category
                double simesh,  // simesfactor at h(alpha)
                double alpha,   // alpha itself
                int    m)       // size of the problem
{
  if (p==0 || simesh==0)
    return 1;
  else
    if (alpha==0)
      return m+1;
    else
    {
      double cat = (simesh/alpha)*p;
      return static_cast<int> (std::ceil(cat));
    }
}

// Calculate the lower bound to the number of false hypotheses
// Implement the algorithm based on the disjoint set structure
std::vector<int> findDiscoveries(matlab::data::TypedArray<int32_t>& ixp,     // subset indices (starts from 1)
                                 matlab::data::TypedArray<double>&  allp,    // all p-values (no need to be sorted)
                                 double                             simesh,  // simesfactor at h(alpha)
                                 int                                h,       // h(alpha)
                                 int                                z,       // size of the concentration set
                                 double                             alpha,   // alpha
                                 int                                k,       // subset size
                                 int                                m)       // size of the problem
{
  // calculate categories for the p-values
  std::vector<int> cats(k);
  for (int i=0; i<k; i++)
  {
    cats[i] = getCategory(allp[ixp[i]-1], simesh, alpha, m);
  }
  
  // find the maximum category needed
  int maxcat  = std::min(z-m+h+1, k);
  int maxcatI = 0;
  for (int i=k-1; i>=0; i--)
  {
    if (cats[i] > maxcatI)
    {
      maxcatI = cats[i];
      if (maxcatI >= maxcat) break;
    }
  }
  maxcat = std::min(maxcat, maxcatI);
  
  // prepare disjoint set data structure
  std::vector<int> parent(maxcat+1);
  std::vector<int> lowest(maxcat+1);
  std::vector<int> rank(maxcat+1,0);
  for (int i=0; i<=maxcat; i++)
  {
    parent[i] = i;
    lowest[i] = i;
  }

  // The algorithm proper. See pseudocode in paper
  std::vector<int> discoveries(k+1,0);
  int lowestInPi;
  for (int i=0; i<k; i++)
  {
    if (cats[i] <= maxcat)
    {
      lowestInPi = lowest[Find(cats[i], parent)];
      if (lowestInPi == 1)
      {
        discoveries[i+1] = discoveries[i]+1;
      }
      else
      {
        discoveries[i+1] = discoveries[i];
        Union(lowestInPi-1, Find(cats[i], parent), parent, lowest, rank);
      }
    }
    else
    {
        discoveries[i+1] = discoveries[i];
    }
  }
  
  return discoveries;
}
