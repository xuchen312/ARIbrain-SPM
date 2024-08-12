/* 
 * spm_aribrain_tdp
 * clusterlist = spm_aribrain_tdp(m,dims,maskI,indexp,ordp,rankp,conn,p);
 * Compute maximal clusters using adaptive thresholding algorithm.
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "spm_aribrain.h"

#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <cmath>

class MexFunction : public matlab::mex::Function {
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
        
        
        
        
        // specify inputs
        int       m = inputs[0];
        int    dims = inputs[1];
        int   maskI = inputs[2];
        int  indexp = inputs[3];
        int    ordp = inputs[4];
        int   rankp = inputs[5];
        double conn = inputs[6];
        double    p = inputs[7];  //std::vector<double> p_vec(p, p + sizeof(p)/sizeof(p[0]));
        
        
        
        
        
        
        
        // h & simes(h)
        int         halpha = findHalpha(jumpalpha, alpha, m);
        double simeshalpha = simesfactor[halpha+1];
        
        // find the adjacency list
        std::vector<std::vector<int>> adj = findAdjList(maskI, as.integer(indexp-1), dims, m, conn);
        
        // find all STCs
        auto reslist = findClusters(m, adj, ordp, rankp);
        // auto [CHILD, ROOT, SIZE] = findClusters(m, adj, ordp, rankp);
        
        // estimate TDP bounds
        std::vector<double> tdps = forestTDP(m, halpha, level, simeshalpha, p, ordp, reslist.SIZE, reslist.ROOT, reslist.CHILD);
        
        // make proper preparations
        std::vector<int> stcs = queryPreparation(m, reslist.ROOT, tdps, reslist.CHILD);
        
        // form maximal clusters based on a TDP threshold
        std::list<std::vector<int>> clusterlist = answerQuery(gamma, stcs, ordp, reslist.SIZE, marks, tdps, reslist.CHILD);
        
        // return outputs (maximal clusters): #{outputs} = #{clusters}
        int i = 0;
        for(std::list<std::vector<int>>::iterator it = clusterlist.begin(); it != clusterlist.end(); ++it)
        {
            for(int j = 0; j < (*it).size(); j++)
            {
                outputs[i][j] = (*it)[j];
            }
            i++;
        }
        
        
        
    }

    
    //------------------------- (1) FIND ALL STCS (USING SORTING RANKS) -------------------------//

    // Union function of a disjoint-set data structure based on the "union by size" technique.
    // The disjoint sets will represent the components of a forest, and the data structure is
    // augmented to keep track of the forest root of each component. UnionBySize(i,j) merges
    // the sets S_i and S_j (where S_k denotes the set containing k), and assigns the forestroot
    // of S_i to be the forestroot of the resulting union.
    void UnionBySize(int               i,
                     int               j,
                     std::vector<int>& PARENT,
                     std::vector<int>& FORESTROOT,
                     std::vector<int>& SIZE)
    {
        int irep = Find(i, PARENT);
        int jrep = Find(j, PARENT);
        
        // if i and j are already in the same set
        if (irep == jrep) return;
        
        // i and j are not in same set, so we merge
        int iroot = FORESTROOT[irep];
        int jroot = FORESTROOT[jrep];
        if (SIZE[iroot] < SIZE[jroot])
        {
            PARENT[irep] = jrep;
            FORESTROOT[jrep] = iroot;
        }
        else
        {
            PARENT[jrep] = irep;
        }
        SIZE[iroot] += SIZE[jroot];
    }
    
    // Compute all supra-threshold clusters (STCs)
    struct clusters
    {
        std::list<std::vector<int>> CHILD;
        std::vector<int> ROOT;
        std::vector<int> SIZE;
    };
    
    auto findClusters(int                            m,     // number of nodes
                      std::vector<std::vector<int>>& ADJ,   // a list of neighbours for all nodes (unsorted!)
                      std::vector<int>&              ORD,   // sorted orders for non-decreasing p-values
                      std::vector<int>&              RANK)  // sorting ranks for all p-values
    {
        // initialize output (1): a list of children for all nodes
        std::list<std::vector<int>> CHILD(m);
        // initialize output (2): a list of sizes of subtrees
        std::vector<int> SIZE(m, 1);
        // initialize output (3): a list of forest roots
        std::list<int> ROOT;
        
        // initialize a child list for a single node
        std::list<int> CHD;
        // // using std::forward_list instead of std::list (Note: must include header <forward_list>)
        // std::forward_list<int> CHD;
        
        // prepare disjoint set data structure
        std::vector<int> PARENT, FORESTROOT;
        PARENT.reserve(m);
        FORESTROOT.reserve(m);
        for (int i = 0; i < m; i++)
        {
            PARENT.push_back(i);
            FORESTROOT.push_back(i);
        }
        
        // loop through all nodes in the ascending order of p-values
        // // for C++11 and 0-based ORD, the use of v is not needed
        // for(int i : ORD)
        for (int i = 0; i < m; i++)
        {
            int v = ORD[i]-1;
            
            // find neighbours for node with the ith smallest p-value
            std::vector<int> IDS = ADJ[v];
            
            // loop through all its neighbours
            for (int j = 0; j < IDS.size(); j++)
            {
                if (RANK[IDS[j]-1] < i+1)  // check if the neighbour has a smaller rank
                {
                    int jrep = Find(IDS[j]-1, PARENT);  // representative of the tree
                    int    w = FORESTROOT[jrep];        // forest root of the tree
                    
                    if (v != w)
                    {
                        // Merge S_v and S_w=S_{jrep}
                        UnionBySize(v, jrep, PARENT, FORESTROOT, SIZE);
                        
                        // put a heavy child in front (using std::list)
                        if (CHD.empty() || SIZE[CHD.front()] >= SIZE[w])
                        {
                            CHD.push_back(w);
                            // CHD.insert_after(CHD.begin(), w);  // for std::forward_list
                        }
                        else
                        {
                            CHD.push_front(w);
                            // CHD.push_back(CHD[0]);  // for std::vector
                            // CHD[0] = w;
                        }
                    }
                }
            }
            
            // update child list
            CHILD[v] = std::vector<int>(CHD.begin(), CHD.end());
            CHD.resize(0);
        }
        
        // find forest roots
        for (int i = 0; i < m; i++)
        {
            if (PARENT[i] == i)
            {
                ROOT.push_back(FORESTROOT[i]);
            }
        }
        
        return clusters{ CHILD, std::vector<int>(ROOT.begin(), ROOT.end()), SIZE };
    }


    //------------------------- (2) COMPUTE TDPS FOR ALL STCS -------------------------//

    // Iterative post-order traversal to find descendants of v (including v).
    // Note: When we pop a vertex from the stack, we push that vertex again as a value and
    // then all its children in reverse order on the stack. If we pop a value, it means that
    // all its children have been fully explored and added to the descendants, so we append
    // the current value to the descendants too.
    std::vector<int> descendants(int                          v,      // node v (0:m-1)
                                 std::vector<int>&            SIZE,   // subtree sizes for all nodes
                                 std::list<std::vector<int>>& CHILD)  // a list of children for all nodes
    {
        std::vector<int> DESC(SIZE[v], 0);
        
        int len = 0;          // track the number of found descendants
        int top = SIZE[v]-1;  // track the top of the stack
        DESC[top] = v;        // push v on the stack (stack grows from right to left)
        
        while (top < DESC.size())  // while stack is non-empty
        {
            v = DESC[top];       // pop the top of the stack
            top++;
            if (v < 0)           // if ~v is a value, append ~v to the descendants
            {
                DESC[len] = ~v;  // bitwise negation is used as minus doesn't work for zero
                len++;
            }
            else
            {
                top--;
                DESC[top] = ~v;  // push v as a value
                
                // push all children in reverse order
                std::vector<int> CHD = *(CHILD.begin()+v);
                for (int j = CHD.size() - 1; j >= 0; j--)
                {
                    top--;
                    DESC[top] = CHD[j];
                }
            }
        }
        
        return DESC;
    }

    // Compute the TDP bounds of the heavy path starting at v
    void heavyPathTDP(int                          v,       // start of the heavy path (0:m-1)
                      int                          par,     // parent of v (-1 indicates no parent)
                      int                          m,       // number of all nodes
                      int                          h,       // h(alpha)
                      double                       alpha,   // alpha
                      double                       simesh,  // simesfactor at h(alpha)
                      std::vector<double>&         P,       // all p-values (unsorted!)
                      std::vector<int>&            SIZE,    // subtree sizes for all nodes
                      std::list<std::vector<int>>& CHILD,   // a list of children for all nodes
                      std::vector<double>&         TDP)     // TDP bounds
    {
        // Rcpp::IntegerVector HP = descendants(v, SIZE, CHILD);
        // for (int i = 0; i < HP.size(); i++)
        // {
        //     HP[i]++;
        // }
        std::vector<int> HP = descendants(v, SIZE, CHILD) + 1;
        std::vector<int> NUM = findDiscoveries(HP, P, simesh, h, alpha, HP.size(), m);
        
        while (true)  // walk down the heavy path
        {
            // check if v represents an STC
            if (par == -1 || P[v] != P[par])
            {
                TDP[v] = ((double) NUM[SIZE[v]]) / ((double) SIZE[v]);
            }
            else
            {
                TDP[v] = -1;  // invalid STCs get TDP of -1
            }

            // check if v is a leaf
            if (SIZE[v] == 1) break;
            
            // update v & its parent
            par = v;
            std::vector<int> CHD = *(CHILD.begin()+v);
            v = CHD[0];
        }
    }

    // Find the start of every heavy path & compute the TDPs of that heavy path
    // start of heavy path: 1) root of F;
    //                      2) non-root node that is not the 1st heavy child
    std::vector<double> forestTDP(int                          m,       // number of all nodes
                                  int                          h,       // h(alpha)
                                  double                       alpha,   // alpha
                                  double                       simesh,  // simesfactor at h(alpha)
                                  std::vector<double>&         P,       // all p-values (unsorted!)
                                  std::vector<int>&            SIZE,    // subtree size for all nodes
                                  std::vector<int>&            ROOT,    // all roots of the forest
                                  std::list<std::vector<int>>& CHILD)   // a child list for all nodes
    {
        std::vector<double> TDP(m);
        
        // loop through all roots
        for (int i = 0; i < ROOT.size(); i++)
        {
            heavyPathTDP(ROOT[i], -1, m, h, alpha, simesh, P, SIZE, CHILD, TDP);
        }
        // loop through all nodes
        for (int i = 0; i < m; i++)
        {
            std::vector<int> CHD = *(CHILD.begin()+i);
            for (int j = 1; j < CHD.size(); j++)
            {
                heavyPathTDP(CHD[j], i, m, h, alpha, simesh, P, SIZE, CHILD, TDP);
            }
        }

        return TDP;
    }


    //------------------------- (3) PREPARE ADMISSIBLE STCS -------------------------//

    // construct a comparator for the below sorting step
    struct compareBy
    {
        std::vector<double>& value;
        compareBy(std::vector<double>& val) : value(val) {}
        bool operator() (int i, int j) {return value[i] < value[j];}
    };

    // Set up ADMSTC: a list of representative of admissible STCs
    std::vector<int> queryPreparation(int                          m,      // number of vertices
                                      std::vector<int>&            ROOT,   // all roots of the forest
                                      std::vector<double>&         TDP,    // all TDP bounds
                                      std::list<std::vector<int>>& CHILD)  // a children list for all vertices
    {
        std::vector<int> ADMSTC;  // a vector of representatives of admissible STCs
        ADMSTC.reserve(m);
        std::vector<double> STACK;
        STACK.reserve(m*2);
        
        // loop through all roots
        for (int i = 0; i < ROOT.size(); i++)
        {
            STACK.push_back(ROOT[i]);  // walk down the forest from ROOT[i]
            STACK.push_back(-1);       // maximum seen TDP on the path to ROOT[i]: non-existent
            while (STACK.size() > 0)
            {
                double q = STACK.back();  // maximum seen TDP on the path to v
                STACK.pop_back();
                int    v = int(STACK.back());
                STACK.pop_back();
                
                // check if v has higher TDP than its ancestors
                if (TDP[v] > q) ADMSTC.push_back(v);  // note: q>=-1 & invalid STCs have TDP=-1
                
                std::vector<int> CHD = *(CHILD.begin()+v);
                for (int j = 0; j < CHD.size(); j++)
                {
                    STACK.push_back(CHD[j]);
                    STACK.push_back(std::max(TDP[v], q));
                }
            }
        }
        
        // sort ADMSTC in ascending order of TDP using the comparator
        std::sort(ADMSTC.begin(), ADMSTC.end(), compareBy(TDP));

        return ADMSTC;
    }


    //-------------------------- (4) FORM CLUSTERS USING gamma --------------------------//

    // Find leftmost index i in ADMSTC such that TDP[ADMSTC[i]] >= g
    // return size(ADMSTC) if no such index exists;
    // run linear search & binary search in parallel;
    // gamma>=0 is needed because inadmissible STCs have been assigned TDP -1.
    int findLeft(double               gamma,    // a TDP threshold (gamma)
                 std::vector<int>&    ADMSTC,   // a list of all admissible vertices (sorted on TDP)
                 std::vector<double>& TDP)      // all TDP bounds
    {
        int right = ADMSTC.size();
        int   low = 0;
        int  high = right;
        while (low < high)
        {
            int mid = (low+high)/2;  // (1) binary search part (using integer division)
            if (TDP[ADMSTC[mid]] >= gamma)
            {
                high = mid;
            }
            else
            {
                low = mid + 1;
            }
            
            right--;                 // (2) linear search part
            // no need to guard against right<0 as right>=0 will always be true
            if (TDP[ADMSTC[right]] < gamma) return (right+1);
        }
        
        return low;
    }

    // Answer the query, i.e., find maximal STCs under the TDP condition.
    // gamma>=0 is needed because inadmissible STCs have been assigned TDP -1.
    std::list<std::vector<int>> answerQuery(double                       gamma,
                                            std::vector<int>&            ADMSTC,
                                            std::vector<int>&            SIZE,
                                            std::vector<int>&            MARK,
                                            std::vector<double>&         TDP,
                                            std::list<std::vector<int>>& CHILD)
    {
        if (gamma < 0) gamma = 0;  // constrain TDP threshold gamma to be non-negative
        
        // initialise output: a list of sorting rank vectors for all clusters
        std::list<std::vector<int>> ANS;
        
        int left = findLeft(gamma, ADMSTC, TDP);
        
        for (int i = left; i < ADMSTC.size(); i++)
        {
            if (MARK[ADMSTC[i]] == 0)
            {
                // append a cluster to ANS
                std::vector<int> DESC = descendants(ADMSTC[i], SIZE, CHILD);
                ANS.push_back(DESC);
                // mark the corresponding voxels
                for (int j = 0; j < DESC.size(); j++)
                {
                    MARK[DESC[j]] = 1;
                }
            }
        }
        
        // clear marks back to 0
        for(std::list<std::vector<int>>::iterator it = ANS.begin(); it != ANS.end(); ++it)
        {
            for(int j = 0; j < (*it).size(); j++)
            {
                MARK[(*it)[j]] = 0;
            }
        }
        
        return ANS;
    }

    // Counting sort in descending order of cluster sizes.
    std::vector<int> counting_sort(int               n,          // #{clusters}
                                   int               maxid,      // max(cluster size)
                                   std::vector<int>& CLSTRSIZE)  // unsorted cluster sizes
    {
        // initialise output sorted indices for descending cluster sizes
        std::vector<int> SORTED(n, 0);
        std::vector<int> COUNT(maxid+1, 0);
        
        // store count of each cluster size
        for (int i = 0; i < n; i++)
        {
            COUNT[CLSTRSIZE[i]]++;
        }
        
        // find cumulative frequency
        for (int i = maxid; i > 0; i--)
        {
            COUNT[i-1] += COUNT[i];
        }
        
        for (int i = 0; i < n; i++)
        {
            SORTED[COUNT[CLSTRSIZE[i]] - 1] = i;
            COUNT[CLSTRSIZE[i]]--;
        }
        
        return SORTED;
    }
    
    
};
