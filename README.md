# ARIbrain-SPM

**ARIbrain-SPM** is an SPM extension for conducting the true discovery proportion (TDP) inference using All-Resolutions Inference (ARI) in neuroimaging.

## Introduction

The free and open-source [Statistical Parametric Mapping (SPM)](https://www.fil.ion.ucl.ac.uk/spm/) software was designed for the analysis of brain imaging data, and the tools for neuroimaging data analyses that are based on SPM are SPM extensions.

Here we propose a new SPM extension **ARIbrain-SPM** that can be used to implement the TDP inference based on the ARI approach.

* The original [ARI](https://doi.org/10.1016/j.neuroimage.2018.07.060) approach was introduced to deal with the spatial specificity paradox. For cluster inference on voxels, there exists at least one activated voxel in a significant cluster, but the acitivation location and amount are unknown. Using ARI, the TDP lower bound, which indicates the activation proportion within each cluster, can be calculated for any arbitrary clusters. The spatial specificity paradox is alliviated when clusters with large enough TDP are identified.

* The [adaptive thresholding algorithm](https://arxiv.org/abs/2206.13587) was developed to further resolve the problem related with the spatial specificity paradox. The ARI approach can be used to estimate the TDP bounds, however, it does not guide how to find the clusters with high enough TDP. Therefore, with a sufficient TDP threshold, the adaptive thresholding algorithm can be applied to find all maximal supra-threshold clusters, for which ARI gives a TDP bound that is at least the chosen threshold. For a typical fMRI dataset, running such an efficient algorithm, with a linearithmic time complexity, simply takes a few seconds.

## References

Rosenblatt, J.D., Finos, L., Weeda, W.D., Solari, A. and Goeman, J.J. (2018). All-Resolutions Inference for brain imaging. *Neuroimage*, 181:786-796. [[Paper](https://doi.org/10.1016/j.neuroimage.2018.07.060)]

Chen, X., Goeman, J.J., Krebs, T.J.P., Meijer, R.J. and Weeda, W.D. (2023). Adaptive Cluster Thresholding with Spatial Activation Guarantees Using All-resolutions Inference. *arXiv*. [[Paper](https://arxiv.org/abs/2206.13587)]

## Any queries?

Please email xuchen312@gmail.com.

