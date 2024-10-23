/* 
 * spm_aribrain_cluster.cpp - Compute the TDP lower bounds for supra-threshold clusters using ARI.
 *
 * The calling syntax is:
 *
 *   tdn = spm_aribrain_cluster(m,k,alpha,simes,ixp,allp,sortp);
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <vector>

#include "spm_aribrain_cluster.h"

class MexFunction : public matlab::mex::Function {
    
    matlab::data::ArrayFactory factory;
    
public:
    
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        // define input arguments (inputs & outputs type - matlab::data::Array)
        int        m = inputs[0][0];  // number of all p-values
        int        k = inputs[1][0];  // number of subset p-values
        double alpha = inputs[2][0];  // significance level
        bool   simes = inputs[3][0];  // (TRUE): Simes' test or (FALSE): Hommel's robust test
        matlab::data::TypedArray<int32_t>  ixp = std::move(inputs[4]);  // subset indices of original p-values (1:m)
        matlab::data::TypedArray<double>  allp = std::move(inputs[5]);  // all original p-values
        matlab::data::TypedArray<double> sortp = std::move(inputs[6]);  // all sorted p-values
        
        // find TDN
        std::vector<double> simesfactors = findsimesfactor(simes, m);
        std::vector<double>    jumpalpha = findalpha(sortp, m, simesfactors, simes);
        int                            h = findHalpha(jumpalpha, alpha, m);
        int                            z = findConcentration(sortp, simesfactors[h], h, alpha, m);
        std::vector<int>     discoveries = findDiscoveries(ixp, allp, simesfactors[h], h, z, alpha, k, m);
        
        // return TDN
        outputs[0] = factory.createScalar(discoveries[k]);
    }
    
};
