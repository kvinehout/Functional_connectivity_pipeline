function [ ts ] = Calculate_Global_Connectivity(PATH_NAME,subject_name,ica_file,good_compi,file_length,lengthcomp,ROIGM,task_name)
%This gets the ORI filtered time series for each cluster componet area
%THIS IS FOR GLOBAL CONNECTIVITY 
%inputs
%path_name: this is the file leading up to the subject folders

%subject_name in ORDER of GLM INPUT FILE WITH LOCAION 
%if multiple conditions then modify subject_name to account for this
%{/con1sub1,/con2sub2}
%good_compi is the list of comps that are used from ICA
%file_length is the lenght of FMRI data collected
%ROIGM is varibale to determine if whole ROI or if only GM ROI is used,
%0=whole ROI, 1=ONLY GM

%outputs: 
%TS is the global connectiity timesieres structure


%% main script 
cd(PATH_NAME)
mkdir global_connectivity
%lengthcomp=length(good_compi);
length_subject_name=length(subject_name);
%pre allocate variables 
alltime=zeros(lengthcomp,(file_length*length_subject_name));    %(lengthtrail_P*file_length_T*length_subject_name_P)+(2*lengthtrail_T*file_length_T*(length_subject_name_T-1))+(file_length_T*lengthtrail_T));
% get the time course for each subject for each componet 
count_compi=1;
for subji=1:length_subject_name%compi=good_compi %for each good_compi values
   count=1;%for naming of timeseries
   %begining of each subject set to zero
   all_subject=zeros(1,(length_subject_name*file_length)); 
   for compi=good_compi%subji=1:length_subject_name  
        %get into the directory
       % all_trails=[];
       %filename=sprintf('%s.feat/sub1.results/errts.sub1.nii',subject_name{subji}); %this is the locaion of the subject ICA_AROMA FOLDER
       if ROIGM==0
            timecourse_data=sprintf('fslmeants -i %s -o global_connectivity/%d_timecourse_comp%d_subji_%d_%d_%snewtwo.txt -m %s/mask_comp_whole%d.nii.gz',subject_name{subji},count,compi,subji,lengthcomp,task_name,ica_file,compi);
       elseif ROIGM==1
            timecourse_data=sprintf('fslmeants -i %s -o global_connectivity/%d_timecourse_comp%d_subji_%d_%d_%snewtwo.txt -m %s/mask_comp_whole_GM_%d.nii.gz',subject_name{subji},count,compi,subji,lengthcomp,task_name,ica_file,compi);
       else
           warning('ROIGM not defined correctly')
       end
       unix(timecourse_data);%this line takes up a lot of time --> anyway to speed up? 
       textname=sprintf('global_connectivity/%d_timecourse_comp%d_subji_%d_%d_%snewtwo.txt',count,compi,subji,lengthcomp,task_name);
       fileID=fopen(textname);
       out_str=textscan(fileID,'%f');
       fclose(fileID);
       timecourse=out_str{:,1};
       timecourse=timecourse';
       %make sure file is the correct size
       if length(timecourse)~=file_length
           warnname=sprintf('%s size %d',subject_name{subji},length(timecourse)); 
           warning(warnname) 
       end 
       %get signal of file_length...this step should do nothing
       timecourse=timecourse(1:file_length);
       %demean the signal
       timecourse=timecourse-mean(timecourse);
       max_time=max(timecourse);
       min_time=min(timecourse);
       if max_time==0 && min_time==0
          warnname=sprintf('%d comp, %d subject,gives zero timecourse max value of %d',compi,subji,max_time); 
          warning(warnname) 
       % return
       end
       all_subject(((subji-1)*length(timecourse)+1):(subji*length(timecourse))) = timecourse;
       clear timecourse
       count=count+1;
        % CHECK IS SUBJ IS ZEROS 
        max_time_S=max(all_subject);
        min_time_S=min(all_subject);
        if max_time_S==0 && min_time_S==0
           warnname=sprintf('%s subject zeros',subject_name{subji+1}); 
           warning(warnname) 
        end        
   end   
    alltime(count_compi,1:length(all_subject))=all_subject;
    count_compi=count_compi+1;
    clear all_subject
    compi
end
ts.ts=alltime;
ts.Nsubjects=length_subject_name;
ts.Nnodes=lengthcomp;
ts.NnodesOrig=lengthcomp;
ts.NtimepointsPerSubject=file_length;
ts.TR=2;
ts.Ntimepoints=length(ts.ts);
ts.ts=ts.ts';
save(sprintf('ts_values_global_%s',task_name),'ts')
   
end

