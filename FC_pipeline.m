function [local_conn,global_conn,p_corrected1Z_local,p_corrected1Z_global] = FC_pipeline(ana_DICOM,fmri_DICOM,path_name,Path_ROI,reg_file,BLOCKEVENTALL,group_all,time_file_all,TR,TE,task_name_all,ref,ROIGM,save_name_all,removeTR,skipBET)
%This function calculates the FC pipeline based on DICOM images from the
%MRI
%this can calculate rest or task-based FC, for task-based a task-design
%signal is regressed from the signal
%this calculates local and global FC
%option visualziation of global and local FC
%klvinehout@gmail.com

%% INPUTS AND OUTPUTS TO RUN WHOLE PROGRAM 
% EXTRNAL FILES NEEDED: 
%regression file: task design .txt file 
%list of ROIS 

%INPUTS: (example)
%path_name='/mnt/data/Kaleb_Data/BSFCSC/'; %name to data subject folder
%reg_file='/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/fsl/standard/MNI152_T1_2mm_brain';
%BLOCKEVENTALL=[2,2,2,2,3,0,3,3];% zero for rest, 2 for block, 3 for event 
%time_file_all{1}='/mnt/data/Kaleb_Data/BSFCSC/timefiles/discreate_pedaling_glm.txt';%regress out event design 
%TR=2.000000; %repition time this is in seconds 
%TE=35.0000; %echo time this is in seconds 
%task_name_all{1}='mem'; 
%ref=sprintf('%s.nii.gz',reg_file);
%addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/AIM2')
%ROIGM=1; %this is 1 to use only GM and zero for whole ROI
%save_name_all='Subject1';
%removeTR=4; %remove the first 4 TR
%skipBET=0; %set to 1 if you want to skip BET
%fmri_DICOM='S23_ped/4_pedaling_task_con'; %path after path_name to fMRI DICOM files
%ana_DICOM='S23_ped/2_Ax_FSPGR_3D';   %path after path_name to anatomical DICOM files
%Path_ROI='/Users/kaleb/Documents/DATA/standard_ICA_BM70/melodic_2mm';path to masked images of ROI in standard sapce 
%group_all{1}=[1,2,3,4,5,6]; %used in stats group 1 vs group 2 comparision (subject order of ana_DICOM in group)
%group_all{2}=[7,8,9,10,11,12]; %used in stats group 1 vs group 2 comparision (subject order of ana_DICOM in group)
    
%OUTPUTS: 
%local connectiviy
%local_conn: (size number ROI , number of MRI files)
%global connectiivty 
%global_conn: (size number ROI , number of ROI, number of MRI files  


%% add all matlab files to path 
%This are paths to suplimentary files. 
addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/NIfTI_20140122') %path to nifti files 
addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/dicm2nii') %path to dcm2nii files 
addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/FSLNets') %path to fslnets
ROBEX_PATH='/mnt/data/inerlUsers/ROBEX/build';%path to robex
linda_path='/mnt/data/inerlUsers/LINDA/'; %path to lind
temp_path='/mnt/data/Kaleb_Data/BSFCSC/timefiles/'; %path to timming files  
ica_file='/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/standard_ICA_BM70/melodic_2mm'; %this is the path to the ROIs
glm_path='/Users/kaleb/Documents/MU_Grad/dissertation/AIM2/glm/';
ants_path='/mnt/data/inerlUsers/antsbin/bin/';
%% Sections below use imputs defined above 
%% preposesiing
%Prepocsing to get AFNI files
%only set up for alternativng fmri data collected in +Z alternation 
%NEED each SUBhect/condtion to have OWN folder with OWN anatomical DATA
%for preprocess
ana_file={ana_DICOM};
fmri_file={fmri_DICOM};
[ fin ] = preproc_ana_fmri(ana_file,fmri_file,path_name,removeTR,ROBEX_PATH,skipBET);


%% registration 
ana_file_length=12;
%for registration 
ana_file={sprintf('%s/T1/ana_or_brain',ana_DICOM)};
fMRI_file={sprintf('%s/shortfmri_or_mcft_brain_bias',fmri_DICOM)};
        
[status] = registartion_a2_stroke_control(path_name,reg_file,ana_file,fMRI_file,ana_file_length,linda_path);


%% filtering 
%for rest only run preprcessing steps --> filt fun data
subject_name={sprintf('%s/shortfmri_or_mcft_brain_bias_ICA_AROMA/',fmri_DICOM')};
%so only subjects for that TASK
all_task{1}=[1];   
delvol=0;%removed in preprocess steps
for i=1:length(all_task)
    subject_name_task=subject_name(all_task{i});
    BLOCKEVENT=BLOCKEVENTALL(i); %2 for block desing, 3 for event desgn, value not matter for rest
    for ii=1:length(subject_name_task)
        fMRI_file_name{ii}='denoised_func_data_nonaggr';
    end
    %this is the model file zero and ones, none for rest
    TIMEFILE=time_file_all{i};
    if i==6 %only dont regress for rest
        hemo_nohemo=0;%run for 
        fsfname='restfeat.fsf';
    else
        hemo_nohemo=1;%run for task
        fsfname='taskfeat.fsf';
    end
    subji=i;
    name=task_name_all{i};     
    [fMRI_TIME_all] = feat_hemoR_nohemoR(temp_path,fsfname,path_name,subject_name_task,ref,fMRI_file_name,TR,TE,delvol,TIMEFILE,hemo_nohemo,subji,name,BLOCKEVENT);
    fMRI_TIME_all_task{i}=fMRI_TIME_all;
    clear subject_name_task
    clear fMRI_file_name
end

%% define ROI (optional) 
%defined from ICA file or list of temporially, but not spacially distinct 

%{
%get clusters
ica_file='standard_ICA_BM70';
path_name='/Users/kaleb/Documents/DATA';
min_cluster_size=1000; %THIS IS A 2CM CUBE FOR 2MM FMRI DATA 
standard_IM='/Users/kaleb/Documents/fsl/data/standard/MNI152_T1_1mm_or_brain.nii.gz'

[lengthcomp,IC_comp_cluster] = Define_ROI(ica_file,path_name,min_cluster_size,standard_IM);

%}




%% EXtract connectiivty for all subjects/conditions 
subject_name_file_task={sprintf('%s/shortfmri_or_mcft_brain_bias_ICA_AROMA/funccondata_%s.nii.gz',fmri_DICOM,task_name_all{1})};
%so only subjects for that TASK
all_task{1}=[1];   
lengthcomp=149;
good_compi=[1:lengthcomp]; %use all ROI initally
for taski=1:size(all_task,2)
    subjects=[100:(100+size(all_task{taski},2))]; %list of subject names %list of subject names 100 base is for contorls  200 base is for stroke 
    %get the TR and the file_length from the first subject in task 
    one_subj=subject_name_file_task{all_task{taski}};
    %get the number of frames for FMRI
    hislicen=sprintf('%s/PrintHeader %s/%s | grep Dimens | cut -d '','' -f 4 | cut -d '']'' -f 1',ants_path,path_name,one_subj);
    [s,hislice]=unix(hislicen);
    %get the TR for FMRI data 
    trn=sprintf('%s/PrintHeader %s/%s | grep "Voxel Spac" | cut -d '','' -f 4 | cut -d '']'' -f 1',ants_path,path_name,one_subj);
    [s,tr]=unix(trn);
    %convert character to number
    file_length=str2num(hislice);
    TR=str2num(tr);
    subject_name_file=subject_name_file_task(all_task{taski});
    task_name=task_name_all{taski};
    %gloabl connectivity extraction 
     [ ts ] = Calculate_Global_Connectivity(path_name,subject_name_file,ica_file,good_compi,file_length,lengthcomp,ROIGM,task_name,subjects,TR);
 %save values for all task
    ts_all{taski}=ts;
    %local connectivity extraction 
    %if output here the within TIME data --> can do stats later 
    ppi=0;%run correaltion not PPI
    timmingfile='';%no timming b/c PPI not used
    ppi_subject=zeros(length(subject_name_file)); %array of 1 and zero for each subject in subject_name_file if PPI or correaltion to be used 
    [ subject_mask_avg_new,subject_mask_avg,stdev_sub] = Calculate_Local_Connectivity(path_name,ica_file,subject_name_file,good_compi,file_length,ppi_subject,TR,timmingfile,ROIGM,task_name,subjects);
    %save values for all task
    subject_mask_avg_new_all{taski}=subject_mask_avg_new;
    subject_mask_avg_all{taski}=subject_mask_avg;
    stdev_sub_all{taski}=stdev_sub;
    %clear variables
    clear subject_name_file
    clear ts
    clear subject_mask_avg_new
    clear subject_mask_avg
    clear stdev_sub
    clear ppi_subject
end

%calculate global connectivity 
for i=1:size(all_task,2)
    global_conn(:,:,i)=nets_netmats(ts_all{i},1,'corr'); 
end
%local connectivity 
local_conn=subject_mask_avg_new_all;

%save files
save(sprintf('%s/subject_mask_avg_new_all_%s.mat',path_name,save_name_all),'subject_mask_avg_new_all')
save(sprintf('%s/stdev_sub_all_%s.mat',path_name,save_name_all),'stdev_sub_all')
save(sprintf('%s/ts_all_%s.mat',path_name,save_name_all),'ts_all')

%% visualization 
%number of connectiions to visualize 
good_compi_task=[1:149];
%gloabl connectivity 
     anatomial_file='anatomical_names_cluster_wholeTAL.mat';
    location_file='max_location_whiten_half_brain.mat';
    size_file='total_Cluster_size_whiten_half_brain.mat';
    %RESIZE ROI FOR TASK COMPI
    cd(Path_ROI)
    load(anatomial_file)
    load(location_file)
    load(size_file)
    for ai=1:length(good_compi_task)
        all_net_name_ROI{ai}=all_net_name{good_compi_task(ai)};
    end
    max_value_whole_ROI=max_value_whole(:,good_compi_task);
    total_cluster_size_ROI=total_cluster_size(good_compi_task);
    saveana=sprintf('anatomical_names_cluster_wholeTAL_taskROI_%s.mat',task_name);
    savemax=sprintf('max_location_whiten_half_brain_taskROI_%s.mat',task_name);
    savetotal=sprintf('total_Cluster_size_whiten_half_brain_taskROI_%s.mat',task_name);
    save(saveana,'all_net_name_ROI');
    save(savemax,'max_value_whole_ROI');
    save(savetotal,'total_cluster_size_ROI');
    location_file_ROI=savemax;
    anatomial_file_ROI=saveana;
    size_file_ROI=savetotal;
    [nodes,xyz_max] = visualize_global(location_file_ROI,size_file_ROI,anatomial_file_ROI,global_conn,Path_ROI,save_name_all);
%local connectivity  
     good_compi=[1:149];
     [ total_mask_sum_corr ] = vizualize_within(local_conn,good_compi,Path_ROI,save_name_all);
%% stats
%run GLM matrix
    num_gr1=size(group_all{1},2);
    num_grp2=size(group_all{2},2); 
    GLMname=sprintf('GLM_%s',save_name_all);
    [ matcon ] = make_glm(num_gr1,num_grp2,glm_path,GLMname);
    %get the desgn matrix files from matcon
    design=sprintf('%s.mat',matcon);
    contrast=sprintf('%s.con',matcon);
    grp=sprintf('%s.grp',matcon);
    %run nonparemetic permutation testing 
    %run stats local 
    [p_uncorrected1Z_local,p_corrected1Z_local]=nets_glm_grp(local_conn,design,contrast,grp,0); 
    %run stats global 
    [p_uncorrected1Z_global,p_corrected1Z_global]=nets_glm_grp(global_conn,design,contrast,grp,0); 
end

