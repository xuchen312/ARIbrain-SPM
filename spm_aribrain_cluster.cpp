/* spm_aribrain_cluster
 * c = spm_aribrain_cluster(a,b);
 * Compute the TDP lower bounds for supra-threshold clusters using ARI.
*/
#include "mex.hpp"
#include "mexAdapter.hpp"
#include "spm_aribrain.h"
#include <vector>

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void spm_aribrain_cluster()(ArgumentList outputs, ArgumentList inputs) {
        
        
        
        
        
        
        
        
        
        // h(alpha)
        halpha      = findHalpha(jumpalpha, level, m);
        simeshalpha = simesfactor[halpha+1];
        
        //
        reslist = findClusters(m, adj, ordp, rankp);
        tdps    = forestTDP(m, halpha, level, simeshalpha, p, ordp, reslist$SIZE, reslist$ROOT, reslist$CHILD);
        stcs    = queryPreparation(m, reslist$ROOT, tdps, reslist$CHILD);
        
        //
        clusterlist = answerQuery(gamma, stcs, ordp, reslist$SIZE, marks, tdps, reslist$CHILD);

        
        
        
        
        checkArguments(outputs, inputs);
        const double offSet = inputs[0][0];
        TypedArray<double> doubleArray = std::move(inputs[1]);
        for (auto& elem : doubleArray) {
            elem += offSet;
        }
        outputs[0] = doubleArray;
        
        
        
    }

    
    void checkArguments(ArgumentList outputs, ArgumentList inputs)
    {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        // Get array factory
        ArrayFactory factory;

        // Check offset argument: First input must be scalar double
        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getType() == ArrayType::COMPLEX_DOUBLE ||
            inputs[0].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("First input must be scalar double") }));
        }

        // Check array argument: Second input must be double array
        if (inputs[1].getType() != ArrayType::DOUBLE ||
            inputs[1].getType() == ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Input must be double array") }));
        }
        // Check number of outputs
        if (outputs.size() > 1) {
            matlabPtr->feval(u"error",
                0,
                std::vector<Array>({ factory.createScalar("Only one output is returned") }));
        }
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    // Implementation of Fortune 1989
    std::vector<int> findhull(int                 m,                // length of p
                              std::vector<double> &p)               // p-values (sorted!)
    {
      // intialize output length
      int r;
      std::vector<int> hull(1);
      hull.push_back(1);
      bool notconvex;
      
      // find the hull
      for (int i=2; i<=m; i++)
      {
        if (i==m || (m-1)*(p[i-1]-p[0]) < (i-1)*(p[m-1]-p[0]))
        {
          do {
            r = hull.size()-1;
            if (r>1)
            {
              notconvex = ((i-hull[r-1])*(p[hull[r]-1]-p[hull[r-1]-1]) >= (hull[r]-hull[r-1])*(p[i-1]-p[hull[r-1]-1]));
            } else {
              if (r==1)
              {
                notconvex = (i*p[hull[1]-1] >= hull[1]*p[i-1]);
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

    // Finds the jumps of h(alpha)
    std::vector<double> findalpha(std::vector<double> &p,            // vector of p-values (sorted!)
                                  int                 m,             // length of p
                                  std::vector<double> &simesfactor,  // denominator of local test
                                  bool                simes)         // assume simes yes or no
    {
      // initialize output
      std::vector<double> alpha(m+1);
      
      // initialize
      std::vector<int> hull = findhull(m,p);
      double Dk = 0;
      int    k  = hull.size()-1;
      int    i  = 1;
      
      // algorithm for alpha*
      while (i<=m)
      {
        if (k>1)
          Dk = p[hull[k-1]-1]*(hull[k]-m+i)-p[hull[k]-1]*(hull[k-1]-m+i);
        if (k>1 && Dk<0)
        {
          k--;
        } else {
          alpha[i-1] = simesfactor[i]*p[hull[k]-1]/(hull[k]-m+i);
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

    // Calculates the denominator of the local test
    std::vector<double> findsimesfactor(bool simes,        // assume simes yes or no
                                        int  m)            // number of p-values
    {
      std::vector<double> simesfactor(m+1);  // denominator of simes test
      double multiplier = 0;                 // extra multiplier term needed for non-simes
      simesfactor[0] = 0;
      if (simes)
        for (int i=1; i<=m; i++)
          simesfactor[i] = i;
      else
        for (int i=1; i<=m; i++)
        {
          multiplier += 1.0/i;
          simesfactor[i] = i*multiplier;
        }

      return simesfactor;
    }

    // Calculate the value of h(alpha) for a given alpha
    int findHalpha (std::vector<double> &jumpalpha,   // points where h jumps
                    double              alpha,        // alpha where h is to be evaluated
                    int                 m)            // size of the multiple testing problem
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
    int findConcentration(std::vector<double> &p,           // vector of p-values (sorted!)
                          double              simesfactor,  // simesfactor at h(alpha)
                          int                 h,            // h(alpha)
                          double              alpha,        // alpha itself
                          int                 m)            // size of the problem
      
    {
      // from m-h we increase z until we fulfil the condition
      int z = m-h;
      if (z>0)  // h=m implies z=0
      {
        while ((z<m) & (simesfactor*p[z-1] > (z-m+h+1)*alpha))
        {
          z++;
        }
      }
      return z;
    }

    // Find function for disjoint set data structure
    // (1) Old recursive version
    // int Find(int              x,
    //          std::vector<int> &parent)
    // {
    //   if (parent[x] != x)
    //   {
    //     parent[x] = Find(parent[x], parent);
    //   }
    //
    //   return parent[x];
    // }
    // (2) iterative version (more stable)
    int Find(int              x,
             std::vector<int> &parent)
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
    void Union(int              x,
               int              y,
               std::vector<int> &parent,
               std::vector<int> &lowest,
               std::vector<int> &rank)
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
    int getCategory(double p,               // p-value for which we need the category
                    double simesfactor,     // simesfactor at h(alpha)
                    double alpha,           // alpha itself
                    int    m)               // size of the problem
    {
      if (p==0 || simesfactor==0)
        return 1;
      else
        if (alpha==0)
          return m+1;
        else
        {
          double cat = (simesfactor/alpha)*p;
          return static_cast<int> (std::ceil(cat));
        }
    }

    // Calculates the lower bound to the number of false hypotheses
    // Implements the algorithm based on the disjoint set structure
    std::vector<int> findDiscoveries(std::vector<int>    &idx,         // indices in set I (starts from 1)
                                     std::vector<double> &allp,        // all p-values
                                     double              simesfactor,  // simesfactor at h(alpha)
                                     int                 h,            // h(alpha)
                                     double              alpha,        // alpha
                                     int                 k,            // size of I
                                     int                 m)            // size of the problem
    {
      // calculate categories for the p-values
      std::vector<int> cats(k);
      for (int i=0; i<k; i++)
      {
        cats[i] = getCategory(allp[idx[i]-1], simesfactor, alpha, m);
      }
      
      // find the maximum category needed
      int z       = findConcentration(allp, simesfactor, h, alpha, m);
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
    
    
    
    
    
};
