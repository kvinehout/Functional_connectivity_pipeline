# Functional_connectivity_pipeline
These files are what i use to do the following:

1) Preprocess fMRI data and Anatomical data with:
    a)bias field correction
    b)spatial and temporal filtering [For task-based connectivity the task signal is regressed out]
    c)pre whitening 
    d) ICA-AROMA is used to remove motion artifact from the fMRI data

2)Registration
  a) Anatomical images are registered with ANTs or LINDA (For lesioned brains)
  b) registration transform is applied to fMRI data

3) ROI Selection
  a) the provided ROI are used to obtain 149 regions of interest  
  
  
4) Global Connectivity 
  a) the average time series is extracted for all voxels within each ROI
  b) all the time series can be correlated with full correlation, partial correlation, or global signal regressed correlation
  c) Correlation values are Fisher-Z transformed 

5) Local Connectivity 
  a) All voxels within each ROI are correlated with each other 
  b) The mean of this Correlation is taken to represent Local Connectivity for each ROI
  C) These data are fisher-Z Transformed 


6) Topology Measurements 
	a) Using the Brain Connectivity Toolbox a variety of topology Measurements can be calculated on the global connectivity values 


7) Statistics 
	a) non-parametric permutation tests are performed with FSL randomize software 


