%MAIN script 

%% INPUTS AND OUTPUTS TO RUN WHOLE PROGRAM 
% EXTRNAL FILES NEEDED: 
%regression file: 
%feat filtering file: 
%text file with reggion
%GLM file for ID rest vs task
%GLM file for task SC STATS
%excell with clinical measuments

%INPUTS: 
path_name='/mnt/data/Kaleb_Data/BSFCSC/'; %name to data subject folder
reg_file='/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/fsl/standard/MNI152_T1_2mm_brain';
BLOCKEVENTALL=[2,2,2,2,3,0,3,3];% zero for rest, 2 for block, 3 for event 
time_file_all{1}='/mnt/data/Kaleb_Data/BSFCSC/timefiles/discreate_pedaling_glm.txt';%regress out event design 
TR=2.000000; %this is in seconds 
TE=35.0000;
task_name_all{1}='mem'; 
ref=sprintf('%s.nii.gz',reg_file);
addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/AIM2')
ROIGM=1; %this is 1 to use only GM and zero for whole ROI
save_name_all='Subject1';
removeTR=4; %remove the first 4 TR
skipBET=0; %set to 1 if you want to skip BET
fmri_DICOM='S23_ped/4_pedaling_task_con'; %path after path_name to fMRI DICOM files
ana_DICOM='S23_ped/2_Ax_FSPGR_3D';   %path after path_name to anatomical DICOM files

%OUTPUTS: 
%subject_mask_avg_new_all_%s.mat %saved in path_name folder this is the
%local connectiviy 

%ts_all_%s.mat'  %saved in the path_name folder this is the global
%timeseries 



%% add all matlab files to path 
%for aim 2 files and maltab files 

addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/NIfTI_20140122') %path to nifti files 
addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/dicm2nii') %path to dcm2nii files 
addpath('/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/matlabfiles/Matlab/FSLNets') %path to fslnets
ROBEX_PATH='/mnt/data/inerlUsers/ROBEX/build';%path to robex
linda_path='/mnt/data/inerlUsers/LINDA/'; %path to lind
temp_path='/mnt/data/Kaleb_Data/BSFCSC/timefiles/'; %path to timming files  
ica_file='/mnt/data/Kaleb_Data/HCP_FILES/HCP_cluster_files/standard_ICA_BM70/melodic_2mm'; %this is the path to the ROIs


%% CODE BELOW HERE SHOULD NOT NEED EDIT

%% preposesiing
%Prepocsing to get AFNI files
%only set up for alternativng fmri data collected in +Z alternation 
%NEED each SUBhect/condtion to have OWN folder with OWN anatomical DATA


 %for preprocess
ana_file={ana_DICOM};
fmri_file={fmri_DICOM};

[ fin ] = preproc_ana_fmri(ana_file,fmri_file,path_name,removeTR,ROBEX_PATH,skipBET);

%new anatomical: %s/T1/ana_or_brain
%new fMRI: %s/shortfmri_or_mcft_brain_bias

%% alll lesions on left so NO flipping of brains --> so also should NOT need to flip for LINDA
%ALL LESIONS SHOULD BE ON LEFT SIDE




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


%% ROI DEFN 
%DEFINE FROM STANDARD TASK/REST ICA 

% HERE WANT IS DIFFERENT THAN REST (STROKE VS STROKE) AND (CONTROL VS
% CONTORL) ---> ADD BOTH TO LIST OF INTERESTED ROI 


%IF languae uses PPI: 
    %ID LANAGUE/CONTROL different from rest
    %RUN PPI ON ALL OF THESE AREAS 

%{


%get clusters
ica_file='standard_ICA_BM70';
path_name='/Users/kaleb/Documents/DATA';
min_cluster_size=1000; %THIS IS A 2CM CUBE FOR 2MM FMRI DATA 
standard_IM='/Users/kaleb/Documents/fsl/data/standard/MNI152_T1_1mm_or_brain.nii.gz'


%change Z stat threhold? currently 1.64 (0.05P?) 
%1000 and 1.64 gives 309 ROI 
%1000 andn 1.84 gives 219 
%1000 and 2.54 gives 149 
[lengthcomp,IC_comp_cluster] = cluster_analysis_standard_space_a2(ica_file,path_name,min_cluster_size,standard_IM);


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
    hislicen=sprintf('/mnt/data/inerlUsers/antsbin/bin/PrintHeader %s/%s | grep Dimens | cut -d '','' -f 4 | cut -d '']'' -f 1',path_name,one_subj);
    [s,hislice]=unix(hislicen);
    %get the TR for FMRI data 
    trn=sprintf('/mnt/data/inerlUsers/antsbin/bin/PrintHeader %s/%s | grep "Voxel Spac" | cut -d '','' -f 4 | cut -d '']'' -f 1',path_name,one_subj);
    [s,tr]=unix(trn);
    %convert character to number
    file_length=str2num(hislice);
    TR=str2num(tr);
    subject_name_file=subject_name_file_task(all_task{taski});
    task_name=task_name_all{taski};
    %gloabl connectivity extraction 
     [ ts ] = cluster_node_ori_data_concat_a2_subj(path_name,subject_name_file,ica_file,good_compi,file_length,lengthcomp,ROIGM,task_name,subjects,TR);
 %save values for all task
    ts_all{taski}=ts;
    %local connectivity extraction 
    %if output here the within TIME data --> can do stats later 
    ppi=0;%run correaltion not PPI
    timmingfile='';%no timming b/c PPI not used
    ppi_subject=zeros(length(subject_name_file)); %array of 1 and zero for each subject in subject_name_file if PPI or correaltion to be used 
    [ subject_mask_avg_new,subject_mask_avg,stdev_sub] = correlatin_within_ROI_a2_subj(path_name,ica_file,subject_name_file,good_compi,file_length,ppi_subject,TR,timmingfile,ROIGM,task_name,subjects);
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

%save files
save(sprintf('%s/subject_mask_avg_new_all_%s.mat',path_name,save_name_all),'subject_mask_avg_new_all')
save(sprintf('%s/stdev_sub_all_%s.mat',path_name,save_name_all),'stdev_sub_all')
save(sprintf('%s/ts_all_%s.mat',path_name,save_name_all),'ts_all')

