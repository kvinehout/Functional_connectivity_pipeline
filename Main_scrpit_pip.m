%MAIN script for PIP
%BY Kaleb Vinehout
%Questions: kaleb.vinehout@mu.edu

%% INPUTS AND OUTPUTS TO RUN WHOLE PROGRAM
% EXTRNAL FILES NEEDED: 
%excell with clinical measuments
%dicom DATA
%3rd party matlab FILES
%feat.fsf files
%GLM for stats

%INPUTS: 
%DICOM DATA

%OUTPUTS: 
%p_correctedZ_all stats for global
%p_correctedZ_all stast for local
%BCT topology calc




%% preposesiing
%Prepocsing to get AFNI files
%only set up for alternativng fmri data collected in +Z alternation
%NEED each SUBhect/condtion to have OWN folder with OWN anatomical DATA

%example inputs
path_name='/Users/kaleb/Documents/MU_Grad/dissertation/AIM2/data/C01_raw'; %name to data subject foldre
ana_file={'C01_session1/3_AxFSPGR3D','C01_session2/2_AxFSPGR3D'};     %anatomical for each FMRI file, also name of folder with anatmoical data
fmri_file={'C01_session1/7_Rest','C01_session1/6_Languetask','C01_session1/8_Memorytask','C01_session2/4_pedalingtaskcon','C01_session2/5_pedalingtaskdiscreate'};     %%name of FMRI file for each subject, also name of folder with fmri data %same as subject_name???
removeTR=4; %remove the first 4 TR
ROBEX_PATH='/Users/kaleb/Documents/Matlab/ROBEX';

[ fin ] = preproc_ana_fmri(ana_file,fmri_file,path_name,removeTR,ROBEX_PATH);

%new anatomical: %s/T1/ana_or_brain
%new fMRI: %s/shortfmri_or_mcft_brain_bias

%% alll lesions on left
%ALL LESIONS SHOULD BE ON LEFT SIDE





%% registration
%example inputs
path_name='/Volumes/Samsung_T3/AIM2/data/';
reg_file='/Users/kaleb/Documents/fsl/data/standard/MNI152_T1_2mm_brain';
linda_path='/Users/kaleb/Documents/MU_Grad/dissertation/AIM2/LINDA_v0.2.7/';
ana_file_length=12;
ana_file={'S03_lang/3_AxFSPGR3D/T1/ana_or_brain','S03_ped/3_AxFSPGR3D/T1/ana_or_brain','S03_ped/3_AxFSPGR3D/T1/ana_or_brain','S03_ped/3_AxFSPGR3D/T1/ana_or_brain','S03_ped/3_AxFSPGR3D/T1/ana_or_brain','S04_lang/3_AxFSPGR3D/T1/ana_or_brain','S04_lang/3_AxFSPGR3D/T1/ana_or_brain','S04_lang/3_AxFSPGR3D/T1/ana_or_brain','S04_ped/6_AxFSPGR3D/T1/ana_or_brain','S04_ped/6_AxFSPGR3D/T1/ana_or_brain','S04_ped/6_AxFSPGR3D/T1/ana_or_brain','S05_lang/3_AxFSPGR3D/T1/ana_or_brain','S05_lang/3_AxFSPGR3D/T1/ana_or_brain','S05_lang/3_AxFSPGR3D/T1/ana_or_brain','S05_ped/3_AxFSPGR3D/T1/ana_or_brain','S05_ped/3_AxFSPGR3D/T1/ana_or_brain','S05_ped/3_AxFSPGR3D/T1/ana_or_brain','S05_ped/3_AxFSPGR3D/T1/ana_or_brain','S06_lang/3_AxFSPGR3D/T1/ana_or_brain','S06_lang/3_AxFSPGR3D/T1/ana_or_brain','S06_lang/3_AxFSPGR3D/T1/ana_or_brain','S06_ped/4_AxFSPGR3D/T1/ana_or_brain','S06_ped/4_AxFSPGR3D/T1/ana_or_brain','S06_ped/4_AxFSPGR3D/T1/ana_or_brain','S06_ped/4_AxFSPGR3D/T1/ana_or_brain','S08_lang/3_AxFSPGR3D/T1/ana_or_brain','S08_lang/3_AxFSPGR3D/T1/ana_or_brain','S08_lang/3_AxFSPGR3D/T1/ana_or_brain','S08_ped/4_AxFSPGR3D/T1/ana_or_brain','S08_ped/4_AxFSPGR3D/T1/ana_or_brain'};
fMRI_file={'S03_lang/8_Languagetask/shortfmri_or_mcft_brain_bias','S03_ped/5_pedalingtaskcon/shortfmri_or_mcft_brain_bias','S03_ped/6_pedalingtaskdis/shortfmri_or_mcft_brain_bias','S03_ped/7_pedalingNP/shortfmri_or_mcft_brain_bias','S03_ped/8_pedalingP/shortfmri_or_mcft_brain_bias','S04_lang/4_Memorytask/shortfmri_or_mcft_brain_bias','S04_lang/5_Rest/shortfmri_or_mcft_brain_bias','S04_lang/7_Languagetask/shortfmri_or_mcft_brain_bias','S04_ped/3_pedalingcon/shortfmri_or_mcft_brain_bias','S04_ped/4_pedalingdis/shortfmri_or_mcft_brain_bias','S04_ped/5_pedalingleft/shortfmri_or_mcft_brain_bias','S05_lang/4_Memorytask/shortfmri_or_mcft_brain_bias','S05_lang/5_Rest/shortfmri_or_mcft_brain_bias','S05_lang/7_Languagetask/shortfmri_or_mcft_brain_bias','S05_ped/4_pedalingcon/shortfmri_or_mcft_brain_bias','S05_ped/5_pedalingdis/shortfmri_or_mcft_brain_bias','S05_ped/6_pedalingNP/shortfmri_or_mcft_brain_bias','S05_ped/7_pedalingP/shortfmri_or_mcft_brain_bias','S06_lang/4_Memorytask/shortfmri_or_mcft_brain_bias','S06_lang/7_Languagetask/shortfmri_or_mcft_brain_bias','S06_lang/5_Rest/shortfmri_or_mcft_brain_bias','S06_ped/5_pedalingcon/shortfmri_or_mcft_brain_bias','S06_ped/6_pedalingdis/shortfmri_or_mcft_brain_bias','S06_ped/7_pedalingnp/shortfmri_or_mcft_brain_bias','S06_ped/8_pedalingp/shortfmri_or_mcft_brain_bias','S08_lang/5_Rest/shortfmri_or_mcft_brain_bias','S08_lang/9_Languagetask/shortfmri_or_mcft_brain_bias','S08_lang/10_Memorytask/shortfmri_or_mcft_brain_bias','S08_ped/6_pedalingtaskcon/shortfmri_or_mcft_brain_bias','S08_ped/7_pedalingtaskdis/shortfmri_or_mcft_brain_bias'};                                                                          
[status] = registartion_a2_stroke_control(path_name,reg_file,ana_file,fMRI_file,ana_file_length,linda_path);




%% filtering
    TIMEFILE='/Users/kaleb/Documents/MU_Grad/dissertation/AIM2/glm/memory_glm.txt'; %not used for rest
    TR=2.000000; %this is in seconds
    TE=35.0000;
    delvol=0;
    temp_path='/Users/kaleb/Documents/MU_Grad/dissertation/AIM2/tmp/';
    fsfname='restfeat.fsf'; %could be restfeat.fsf or 'taskfeat.fsf';
    path_name='/Users/kaleb/Documents/MU_Grad/dissertation/';
    ref='/Users/kaleb/Documents/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz';
     subji=1;
    name='dataname';
     BLOCKEVENT=2;%2 for block desing, 3 for event desgn, value not matter for rest
  [fMRI_TIME_all] = feat_hemoR_nohemoR(temp_path,fsfname,path_name,subject_name_task,ref,fMRI_file_name,TR,TE,delvol,TIMEFILE,hemo_nohemo,subji,name,BLOCKEVENT);
    
%% ROI DEFN 
%DEFINE FROM STANDARD TASK/REST ICA 

%THIS STEP IS ALREADY DONE AND ROI ARE INCLUDED
%get clusters
ica_file='standard_ICA_BM70';
path_name='/Users/kaleb/Documents/DATA';
min_cluster_size=1000; %THIS IS A 2CM CUBE FOR 2MM FMRI DATA 
standard_IM='/Users/kaleb/Documents/fsl/data/standard/MNI152_T1_1mm_or_brain.nii.gz'


[lengthcomp,IC_comp_cluster] = cluster_analysis_standard_space_a2(ica_file,path_name,min_cluster_size,standard_IM);



%% EXtract connectiivty for all subjects/conditions 

%global connectiviy extration 
subject_name_file={'C01_session1/7_Rest/shortfmri_or_mcft_brain_bias_ICA_AROMA/funccondata.nii.gz','C01_session1/6_Languetask/shortfmri_or_mcft_brain_bias_ICA_AROMA/funccondata.nii.gz','C01_session1/8_Memorytask/shortfmri_or_mcft_brain_bias_ICA_AROMA/funccondata.nii.gz','C01_session2/4_pedalingtaskcon/shortfmri_or_mcft_brain_bias_ICA_AROMA/funccondata.nii.gz','C01_session2/5_pedalingtaskdiscreate/shortfmri_or_mcft_brain_bias_ICA_AROMA/funccondata.nii.gz','S01_lang/6_Rest/shortfmri_or_mcft_ICA_AROMA/funccondata.nii.gz','S01_lang/9_Memorytask/shortfmri_or_mcft_ICA_AROMA/funccondata.nii.gz','S01_ped/4_pedalingtask1/shortfmri_or_mcft_ICA_AROMA/funccondata.nii.gz','S01_ped/5_pedalingtask2/shortfmri_or_mcft_ICA_AROMA/funccondata.nii.gz'};     %%name of FMRI file for each subject, also name of folder with fmri data %same as subject_name??? 
ica_file='/Users/kaleb/Documents/DATA/standard_ICA_BM70/melodic_2mm'; %path to melodic file
lengthcomp=149;
path_name='/Volumes/Samsung_T3/AIM2/data'; 
ROIGM=1; %this is 1 to use only GM and zero for whole ROI
good_compi=[1:lengthcomp]; %use all ROI initally
task_name='BOTOX TASK';
file_length=190;  %length of MRI file volumes
%gloabl connectivity extraction 
[ ts ] = cluster_node_ori_data_concat_a2(path_name,subject_name_file,ica_file,good_compi,file_length,lengthcomp,ROIGM,task_name);
   

%local connectivity extraction 
good_compi=[1:lengthcomp]; %use all ROI initally
%if output here the within TIME data --> can do stats later
ppi=0;%run correaltion not PPI
ppi_subject= %array of 1 and zero for each subject in subject_name_file if PPI or correaltion to be used
TR=2; %if tr is 2 seconds
ROIGM=1; %this is 1 to use only GM and zero for whole ROI
task_name='BOTOX TASK';
[subject_mask_avg_new,subject_mask_avg,stdev_sub] = correlatin_within_ROI_a2(path_name,ica_file,subject_name_file,good_compi,file_length,ppi_subject,TR,timmingfile,ROIGM,task_name);
     
    %% STATS AFTER TASK ROI DEFINED AND REGIONS CONNECTIVITY VALUES FOUND 
    
    %define stats
    matcon=matcon_all{taski};
    design=sprintf('%s.mat',matcon);
    contrast=sprintf('%s.con',matcon);
    grp=sprintf('%s.grp',matcon);
    subjectused=taskorder{taski};
    
    %RUN STATS ON LOCAL CONNECIVITY
    good_compi=%this is the ROI that you are interested in, or could run on all 149 [ROI EX: good_compi=1:149]
    %limit ts by   good_compi
    subject_mask_avg_new_short=subject_mask_avg_new_task(:,good_compi);
    %local connectivity ON TASK NETWORK
   [p_uncorrected1Z_within,p_corrected1Z_within]=nets_glm_grp(subject_mask_avg_new_short,design,contrast,grp,0);
    %this is the stats of the local connectivity
    
    %Global connectiivty ON TASK NETWORK
    ts_task=ts_all{taski};
    %limit ts by   good_compi
    ts_CS=ts_task;
    ts_CS.ts=ts_task(:,good_compi);
    ts_CS.Nnodes=length(good_compi);
    ts_CS.NnodesOrig=length(good_compi);
    ts_CS.Ntimepoints=length(ts_CS.ts);
    corrset=3; %set for full correaltion values
        [sig_corr_all,sig_uncorr_all,p_uncorrectedZ,p_correctedZ,netmatsZ,all_full] = fslNets_pro_a2(ts_CS,design,contrast,group,grp,corrset,timmingfile,TR,fmri_length);
    %save varaibles for all tasks
        %these are the stats of the global netwroks
        sig_corr_all_all{taski}=sig_corr_all;
    sig_uncorr_all_all{taski}=sig_uncorr_all;
    p_uncorrectedZ_all{taski}= p_uncorrectedZ;
    p_correctedZ_all{taski}=p_correctedZ;
    netmatsZ_all{taski}=netmatsZ;
    all_full_all{taski}=all_full;
    
    
    %GRAPH THEORY METRICS ON WHOLE BRAIN (all nodes)
    %get the correlations values for all connections
        ts_CS_all=ts_task;
    corrset=3; %set for full correaltion values
        [sig_corr_all,sig_uncorr_all,p_uncorrectedZ,p_correctedZ,netmatsZ,all_full] = fslNets_pro_a2(ts_CS_all,design,contrast,group,grp,corrset,timmingfile,TR,fmri_length);
    re_netmatsZ=reshape(netmatsZ,subjectused,lengthcomp,lengthcomp);
    %perform topolgy calculations
    %GRAPH THEORY METRICS  for subjects
        [median_all,SD_all,P,all_dis,stats ] = conection_distance(re_netmatsZ,xyz_max);
    %maybe see if distance between task ROI is on average different
        %depedning on the task?
        threshold=0; %threshold value, min value to count, so zero is all values (positive)
    sig_corr_all_mean=sig_corr_all(1:length(group),:,:);%limit to only mean cases
    [lambda_all,GE_all,LE_all,Clus_all,Spos_all,Sneg_all,vpos_all,vneg_all,cen_all,r_all,cen_all_edge] = BCT_func(sig_corr_all_mean,threshold);
    lambda_all_all{taski}=lambda_all;
    GE_all_all{taski}=GE_all;
    LE_all_all{taski}=LE_all;
    Clus_all_all{taski}=Clus_all;
    Spos_all_all{taski}=Spos_all;
    Sneg_all_all{taski}=Sneg_all;
    vpos_all_all{taski}=vpos_all;
    vneg_all_all{taski}=vneg_all;
    cen_all_all{taski}=cen_all;
    r_all_all{taski}=r_all;
    cen_all_edge_all{taski}=cen_all_edge;
    clear  subject_name_task

    
    
    
    
 %% CLINICAL correlations 

for taski=1:lenght_task 
     clinical_file ='clinicalmeasuresdata'; %change name of clinical file depeding on task
    sig_thresClin=0.05; %set the threshold for comparision P vlaue of 0.05 
    %define all_resutls both local/global/top
    all_resutls=
    %define all_resutls_name both local/global/top
    all_resutls_name=
    
    groupST=
    
    group_name=
    
    subject_nameST=
    
    %run correlations for task ROI local/global for EACH task
    [ R_all_group,P_all_group,RL_all_group,RU_all_group,sig_all_clinic_R] = clinicalmeasures_a2(path_nameCL,subject_name,all_results,clinical_file,results_name,sig_thres,groupST,group_name);

end
 
