function [ subject_mask_avg_new,subject_mask_avg,stdev_sub] = Calculate_Local_Connectivity(path_name,ica_file,subject_name,good_compi,fmri_length,ppi_subject,TR,timmingfile,ROIGM,task_name)
% this finds the local conectiivty foreach subject 
%INPUTS:
%path_name:  this is the file leading up to the subject folders
%ica_file: this is used as the locaion of the MASK ROI defined files
%subject_name: list of subject folder in order defined by stats matrix
%design this is the design file from fsl glm
%contrast this is the contrast file from fsl glm
%grp this is the .grp file made from the stats file
%ROIGM is to use full or GM only ROI, 0=full, 1=only GM 
%file_length: this is the number of volumes of the fMRI data collected
%good_compi: list of the ROI that want to run stats on, based on cluster ROI defined regions
%ppi: this is 1 for ppi correlation, zero for normal correaltion 
%TR=te in secs example: TR=2;
%timmingfile this is the text file with zero/one/-1 extra for the task ONLY
%USED FOR PPI ANALYSIS 


%OUTPUTS: 
 %p_uncorrected1Z_within: uncorretected GROUP STATS on mean data: (number of comparisions , groups)
 %p_corrected1Z_within:  corretected GROUP STATS on mean data: (number of comparisions , groups)
 %subject_mask_avg_new: SAME AS subject_mask_avg
 %subject_mask_avg: for each subject, ROI give mean of all voxels correation values
 %stdev_sub: for each subject, ROI give standard deviation of all voxels correation values


%% main script 
cd(path_name)
mkdir local_connectvity
length_subject_name=length(subject_name);
count_compi=1;
subject_mask_avg=zeros(length_subject_name,length(good_compi));
stdev_sub=zeros(length_subject_name,length(good_compi));
count=1;%for naming of timeseries
for subji=1:length_subject_name  %compi=good_compi 
    % get the data for each subject 
    for compi=good_compi %subji=1:length_subject_name 
                    %define mask name 
                    if ROIGM==0
                        mask_name=sprintf('%s/mask_comp_whole%d.nii.gz',ica_file,compi);
                    elseif ROIGM==1
                        mask_name=sprintf('%s/mask_comp_whole_GM_%d.nii.gz',ica_file,compi);
                    else
                        warning('ROIGM not defined correcly')
                    end
                   %load the 4D time sereis data for CONCAT data 
                % filename=sprintf('%s/%s_ICA_AROMA/%s.feat/sub1.results/errts.sub1.nii',subject_name{subji});
                 subject_name_one=sprintf('%s',subject_name{subji});%,subject_name{subji},filename);
                    %do this for timesereis extraction --> gets around matalb data limits 
                    textfilename=sprintf('local_connectvity/ROI%d_subj_%d_%s.txt',compi,subji,task_name);
                    %this line takes a while
                    unix(sprintf('3dmaskdump -noijk -xyz -mask %s %s > %s',mask_name,subject_name_one,textfilename));
                    %organize this text file into a matlab array for this ROI 
                    fileID = fopen(textfilename);
                    textname=[];
                    unitname='%f';
                    for i=1:(3+fmri_length) %3 offset to ignore voxel locaion data   
                         textname = [textname, ' ', unitname];
                    end
                    C = textscan(fileID,textname);
                    fclose(fileID);
                    %make array (time,voxel)
                    %ignore the first 3 points (x,y,z corr) in each dataset
                    for i=4:(3+fmri_length) 
                        mi=i-3;
                        mask_time_zeros(mi,:)=C{i}; 
                    end
                    %check if all time for one voxel is all zeros--IF
                        %SO THAN IGNORE THIS VOXEL 
                        count2=1;
                    for ii=1:size(mask_time_zeros,2)
                        %if all equal to zero or some other number
                        if max(mask_time_zeros(:,ii))==min(mask_time_zeros(:,ii))
                            %skip this voxel 
                             warning(sprintf('Voxel %d skipped in ROI %d due to constant value of %d',ii,compi,max(mask_time_zeros(:,ii))));
                        else
                            mask_time(:,count2)= mask_time_zeros(:,ii);
                            count2=count2+1;
                        end
                    end
                    %demean the mask signal, based on ALL ROI data
                    mean_mask_time=mask_time-mean(mask_time(:));
                    %get correaltion AVERGAE of ALL voxels in the mask use full correaltion
                    ts_subMask.ts=mean_mask_time;
                    ts_subMask.Nsubjects=1; % is this number of subjects is one?  
                    ts_subMask.Nnodes=size(mask_time,2);
                    ts_subMask.NnodesOrig=size(mask_time,2);
                    ts_subMask.NtimepointsPerSubject=fmri_length;
                    ts_subMask.TR=TR;
                    ts_subMask.Ntimepoints=size(ts_subMask.ts,1);
                    %THIS APPLIES Z TRANSFORM FOR EACH SUBJECT/ROI 
                    %SO ALL subject/ROI are already normalized when put
                    %into larget matrix
                    ppi= ppi_subject(subji);
                    if ppi==0
                        netmats1Z=nets_netmats(ts_subMask,1,'corr'); 
                        netmats1Z(netmats1Z<0)=[];%get rid of correlations less than zero
                        timmingfile=0; %this file is not used
                    elseif ppi==1
                           %this is done once per subject and per ROI 
                         [Onsetfile, Movementfile,RT,dt,fMRI_T0,HParam] = SPM_files_PPI(timmingfile,TR,fmri_length);
                         [PPI_beta,PPI_tfac,PPI_trRV]=wb_gppi(mask_time',Onsetfile, Movementfile, subject_name{subji},RT,dt,fMRI_T0,HParam);
                         netmats1Z=PPI_beta; %not sure if in same order here? 
                    end
                    %apply the fisher Z to this...than run stats this
                    %should be better ---> average and stadnard dev --->
                    %should be contasnt for fisher Z --zscore data 
                    avg_corr_subject=mean(netmats1Z(:));
                    stdev_sub(count,count_compi)=std(avg_corr_subject);
                   subject_mask_avg(count,count_compi)=avg_corr_subject; %needs to be subjecs X elements for sig testing 
                   clear ts_subMask
                   clear netmats1Z
                   clear mask_time
                   clear x
                   clear y
                   clear z
                   clear mask_time_zeros
                   clear mask_time
                   count_compi=count_compi+1;       
    end 
    compi
    count=count+1;
    clear  mask_compi 
    %save compi data
    save_name=sprintf('within_mean_corr_subji_%dZ_%s',subji,task_name);
    save(save_name,'subject_mask_avg') 
end
subject_mask_avg_new=subject_mask_avg;

end




