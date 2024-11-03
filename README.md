# ARIbrain-SPM

**ARIbrain-SPM** is an SPM extension for conducting the true discovery proportion (TDP) inference using All-Resolutions Inference (ARI) in neuroimaging.

## Introduction

The free and open-source [Statistical Parametric Mapping (SPM)](https://www.fil.ion.ucl.ac.uk/spm/) software was designed for the analysis of brain imaging data, and the toolboxes for neuroimaging data analyses that are developed based on SPM are recognised as SPM extensions. Here we propose a new SPM extension **ARIbrain-SPM** that can be used to implement the TDP inference based on the ARI approach.

The original [ARI](https://doi.org/10.1016/j.neuroimage.2018.07.060) approach was introduced to deal with the spatial specificity paradox. That is, for cluster inference on voxels, there exists at least one activated voxel in a significant cluster, but the activation location and amount are unknown. Using ARI, the TDP lower bound, which quantifies the activation proportion within each cluster, can be estimated for any arbitrary cluster. The spatial specificity paradox is alleviated when clusters with large enough TDP are identified.

The [adaptive thresholding algorithm](https://arxiv.org/abs/2206.13587) was developed to further resolve the spatial specificity paradox. The ARI approach can be used to compute the TDP lower bounds for certain clusters, however, it does not guide how to find such clusters with a reasonably high TDP. Therefore, with a sufficient TDP threshold, the adaptive thresholding algorithm can be applied to find all maximal supra-threshold clusters (STCs). Each of these clusters has a TDP equal to or exceeding the given threshold.

Therefore, the ARI-based TDP inference could be separated into two parts:

1. **clusterTDP** - Given any arbitrary clusters, the corresponding TDP lower bounds are estimated using ARI.

2. **tdpCluster** - Given a sufficient TDP threshold, maximal STCs are created using the adaptive thresholding algorithm.

## Installation

### Prerequisites

* Please download and install Matlab. For macOS users, you could edit the ```.bash_profile``` file and add Matlab to the ```PATH``` by appending
  ``` r
  export PATH=/Applications/MATLAB_***.app/bin:$PATH
  ```
  where the installed Matlab version ```MATLAB_***``` could be found by running ```matlabroot``` in Matlab.

* Please download SPM12 and add it to the Matlab search path. You could follow either

    + **HOME -> Set Path -> Add Folder...**
    
    + Run the following line in Matlab
      ``` r
      addpath(genpath('.../spm12'));
      ```
  
### Installing ARIbrain-SPM

* Please download the latest version of ARIbrain-SPM with
  ``` r
  git clone https://github.com/xuchen312/ARIbrain-SPM.git
  ```

* Please add the folder for the ARIbrain-SPM toolbox to the Matlab search path by following either

  + **HOME -> Set Path -> Add Folder...**

  + Run the below script from Matlab console
    ```r
    addpath(genpath('.../ARIbrain-SPM'))
    ```

## Implementation

* For a Linux server user, please connect to the remote server by enaabling X11 forwarding using the following command. This allows users to run graphical applications on the server and display them locally.
  ```r
  ssh -X username@server
  ``` 
  where ```username``` and ```server``` should be replaced with the username and address/hostname of the remote server, respectively.
  
  NOTE: X11 forwarding must also be enabled on the server.

* Navigate to the folder for the ARIbrain-SPM toolbox with
  ```r
  cd .../ARIbrain-SPM
  ```
  
* Launch Matlab, or execute Matlab from the Terminal (command prompt) without the full desktop GUI while still allowing to display graphs with the command
  ```r
  matlab -nodesktop -nosplash
  ```
  
* Conduct the ARIbrain inference by running the function ```spm_aribrain``` with an even number of inputs in the console. For each input variable, the name of it should be first specified, and the value passed in to the function is followed. The function can be created using the following syntax:
  ```r
  spm_aribrain(['xSPM',xSPM,'alpha',alpha,'file',file,'simes',simes,'conn',conn,'tdpth',tdpth])
  ```
  Here, at most 6 input pairs could be specified, where **clusterTDP** and **tdpCluster** require at most 4 and 6 pairs (i.e., 8 and 12 input arguments), respectively.

  - Input variables for **clusterTDP** & **tdpCluster**
  
      + ```xSPM``` an SPM input data structure detailed in ```spm_getSPM.m```. If not specified, it will be computed interactively.
   
      + ```alpha``` a significance level. If not specified, ```alpha = 0.05``` will be used by default.
   
      + ```file``` a character array specifying the file name for saving the result table to an output CSV file. If not specified, the output table will not be saved.
   
      + ```simes``` a logical variable implying the Simes test is conducted. The Hommel robust test is used if set to be ```false```. If not specified, the Simes test will be employed by default.
 
  - Additional parameters for **tdpCluster**
 
      + ```conn``` a connectivity criterion, 6 (face), 18 (edge) or 26 (vertex). If not specified, 18 will be used by default.
   
      + ```tdpth``` a chosen TDP threshold. If not specified, **tdpCluster** will not be performed.
 
  NOTE: **clusterTDP** and **tdpCluster** can not be implemented simultaneously.

* In addition, some outputs of the ARIbrain inference can be returned for interactive exploration of the results in the control panel with, e.g.,
  ```r
  [hReg,xSPM,SPM,TabDat] = spm_aribrain;
  ```

* Alternatively, the above steps could be executed from the Terminal (command prompt) and quit Matlab in the end with, e.g.,
  ```r
  matlab -nodesktop -nosplash -r "cd('.../ARIbrain-SPM'); spm_aribrain; exit"
  ```

## Result Display

The main **ARIbrain-SPM** results are summarised with a result table ```TabDat``` that can be printed on the Matlab console, visualised from the graphics window in SPM, returned to the workspace, and exported to a CSV file. Here, the summary table is highly related to the SPM12 statistics results table, and the summary variable ```TDP``` shows the lower bound of TDP bound, derived using ARIbrain. Examples of such summary tables are as follows:

1. **clusterTDP**
   
    By running ```spm_aribrain```, the standard cluster thresholding embedded in SPM is performed. The TDP lower bounds for all significant clusters are estimated, and the results are returned with a summary table.
   
    ```
    Statistics: p-values adjusted for search volume
    ================================================================================
    set	set	cluster	cluster	cluster	cluster	peak	peak	peak	peak	peak	
    p	c	p(FWE)	p(FDR)	k	TDP	p(FWE)	p(FDR)	T	Z	p(unc)	x,y,z {mm}
    --------------------------------------------------------------------------------
    0.052	18	0.000	0.000	5894	0.302	0.000	0.001	 11.90	 6.30	0.000	 58 -14   4 	
  		                                        0.001	0.002	  9.72	 5.76	0.000	 56 -22  -4 	
  		                                        0.005	0.003	  8.96	 5.54	0.000	 48 -30   0 	
     		0.000	0.000	4039	0.286	0.001	0.002	  9.98	 5.83	0.000	-58 -14   0 	
  		                                        0.002	0.002	  9.42	 5.68	0.000	-54 -30   6 	
  		                                        0.004	0.003	  9.09	 5.58	0.000	-56 -22   4 	
    		0.000	0.000	276	0.000	0.109	0.016	  7.25	 4.96	0.000	 52   2  52 	
  		                                        0.815	0.143	  5.58	 4.24	0.000	 50   4  44 	
    		0.042	0.009	125	0.000	0.277	0.040	  6.60	 4.70	0.000	 18  -4 -14 	
  		                                        1.000	0.782	  4.07	 3.41	0.000	 10  -8 -12 	
    		0.971	0.402	27	0.000	0.783	0.139	  5.65	 4.27	0.000	 10  -2  -2 	
    		0.877	0.307	36	0.000	0.851	0.145	  5.51	 4.21	0.000	-34  12 -24 	
    		0.302	0.061	72	0.000	0.989	0.281	  4.99	 3.94	0.000	-42  28  -2 	
    		0.034	0.009	131	0.000	1.000	0.521	  4.51	 3.67	0.000	-60  16  32 	
  		                                        1.000	0.635	  4.31	 3.56	0.000	-54  14  26 	
  		                                        1.000	0.708	  4.21	 3.49	0.000	-38  12  26 	
    		0.999	0.536	17	0.000	1.000	0.620	  4.36	 3.58	0.000	 40   4 -44 	
    		0.999	0.536	17	0.000	1.000	0.635	  4.30	 3.55	0.000	 -8  -8 -12 	
  		                                        1.000	0.919	  3.68	 3.16	0.001	-16  -8 -14 	
    		0.999	0.536	16	0.000	1.000	0.741	  4.16	 3.46	0.000	-48  -8  46 	
    		1.000	0.667	10	0.000	1.000	0.770	  4.12	 3.44	0.000	 28 -20  -6 	
    		0.994	0.530	21	0.000	1.000	0.770	  4.11	 3.43	0.000	  8  48  36 	
    		0.956	0.401	29	0.000	1.000	0.807	  4.03	 3.38	0.000	  8  12  60 	
  		                                        1.000	0.850	  3.95	 3.34	0.000	  8   0  64 	
  		                                        1.000	0.861	  3.78	 3.22	0.001	 -2   4  62 	
    		1.000	0.667	10	0.000	1.000	0.821	  3.99	 3.36	0.000	-36  12 -40 	
    		1.000	0.667	11	0.000	1.000	0.850	  3.94	 3.33	0.000	 10  20  58 	
    		1.000	0.667	11	0.000	1.000	0.850	  3.93	 3.32	0.000	  6 -32   0 	
    		0.997	0.536	19	0.000	1.000	0.861	  3.80	 3.24	0.001	 10 -12   8 	
    
    table shows 3 local maxima more than 8.0mm apart
    --------------------------------------------------------------------------------
    Height threshold: T = 3.58, p = 0.001 (1.000)
    Extent threshold: k = 10 voxels, p = 0.279 (1.000)
    Expected voxels per cluster, <k> = 9.220
    Expected number of clusters, <c> = 11.71
    FWEp: 7.728, FDRp: 6.533, FWEc: 125, FDRc: 125
    Degrees of freedom = [1.0, 19.0]
    FWHM = 10.2 9.9 9.2 mm mm mm; 5.1 4.9 4.6 {voxels}
    Volume: 2866384 = 358298 voxels = 2905.8 resels
    Voxel size: 2.0 2.0 2.0 mm mm mm; (resel = 116.13 voxels)
    ================================================================================
    ```

2. **tdpCluster**

    For ```tdpth = 0.7```, running ```spm_aribrain('tdpth',0.7)``` would enforce the implementation of adaptive thresholding procedure. The generated maximal clusters and the corresponding TDP bounds are summarised in the result table below.

    ```
    Statistics: p-values adjusted for search volume
    ================================================================================
    set	set	cluster	cluster	cluster	cluster	peak	peak	peak	peak	peak	
    p	c	p(FWE)	p(FDR)	k	TDP	p(FWE)	p(FDR)	T	Z	p(unc)	x,y,z {mm}
    --------------------------------------------------------------------------------
    1.000	2	0.000	0.000	2292	0.700	0.000	0.001	 11.90	 6.30	0.000	 58 -14   4 	
                                                    0.001	0.002	  9.72	 5.76	0.000	 56 -22  -4 	
                                                    0.005	0.003	  8.96	 5.54	0.000	 48 -30   0 	
        		0.000	0.000	1642	0.700	0.001	0.002	  9.98	 5.83	0.000	-58 -14   0 	
                                                    0.002	0.002	  9.42	 5.68	0.000	-54 -30   6 	
                                                    0.004	0.003	  9.09	 5.58	0.000	-56 -22   4 	






    table shows 3 local maxima more than 8.0mm apart
    --------------------------------------------------------------------------------
    Height threshold: T = 3.58, p = 0.001 (1.000)
    Extent threshold: k = 10 voxels, p = 0.279 (1.000)
    Expected voxels per cluster, <k> = 9.220
    Expected number of clusters, <c> = 11.71
    FWEp: 7.728, FDRp: 6.533, FWEc: 125, FDRc: 125
    Degrees of freedom = [1.0, 19.0]
    FWHM = 10.2 9.9 9.2 mm mm mm; 5.1 4.9 4.6 {voxels}
    Volume: 2866384 = 358298 voxels = 2905.8 resels
    Voxel size: 2.0 2.0 2.0 mm mm mm; (resel = 116.13 voxels)
    ================================================================================
    ```

## References

Rosenblatt, J.D., Finos, L., Weeda, W.D., Solari, A. and Goeman, J.J. (2018). All-Resolutions Inference for brain imaging. *NeuroImage*, 181:786-796. [[Paper](https://doi.org/10.1016/j.neuroimage.2018.07.060)]

Chen, X., Goeman, J.J., Krebs, T.J.P., Meijer, R.J. and Weeda, W.D. (2023). Adaptive Cluster Thresholding with Spatial Activation Guarantees Using All-resolutions Inference. *arXiv*. [[Paper](https://arxiv.org/abs/2206.13587)]

## Found bugs, or any questions?

Please email xuchen312@gmail.com.

