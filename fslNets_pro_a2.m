function [sig_corr_all,sig_uncorr_all,p_uncorrectedZ,p_correctedZ,netmatsZ,all_full] = fslNets_pro_a2(ts,design,contrast,group,grp,corrset,timmingfile,TR,fmri_length)
%This runs fsl nets stats 
%this does NOT work for data with multiple trials

%INPUT
%ts this is the subject by connection timeseries data
%design this is the design file from fsl glm
%contrast this is the contrast file from fsl glm
%group: this is an array detailing whitch subjects is ts and the stats
    %matrix belong to what group
    %EX: group{1}=[1,3,5,6,9]; these subjects are from one group in design IN
    %ORDER OF DESIGN MATRIX
    % group{2}=[2,4,7,8,10]; these subjects are from another group IN
    %ORDER OF DESIGN MATRIX
    %group{n} for n number of groups is experimnet
%grp this is the .grp file made from the stats file

%OUTPUT (comparisions,number of connecitons, number of connections) 
%value give is 1-P and ONLY 0.95 are given here
% IN THIS FORM: 
%sig_corr_all_P partial correation corrected ONLY SIGNIFCANT CONNECTIONS,
%sig_uncorr_all_P partial correation uncorrected ONLY SIGNIFCANT CONNECTIONS,
%sig_corr_all_F: full correation correted with global signal regression ONLY SIGNIFCANT CONNECTIONS,
%sig_uncorr_all_F: full correation ucorreted with global signal regression ONLY SIGNIFCANT CONNECTIONS,
%sig_corr_all_F1: full correation correted ONLY SIGNIFCANT CONNECTIONS,
%sig_uncorr_all_F1: full correation ucorreted ONLY SIGNIFCANT CONNECTIONS,

%IN THE FORM:(comparisions,number of connecitons, number of connections) 
%value give is 1-P
%p_uncorrectedp5Z_SM: partial correation uncorrected ALL VALUES
%p_correctedp5Z_SM: partial correation corrected ALL VALUES
%p_uncorrectedf1rZ_SM: full correation uncorreted with global signal regression ALL VALUES
%p_correctedf1rZ_SM: full correation correted with global signal regression ALL VALUES
%p_uncorrectedf1Z_SM: full uncorreation correted ALL VALUES
%p_correctedf1Z_SM: full correation correted ALL VALUES

%IN THE FORM: (subjects,[length(ROI) X length(ROI)] )
%netmats1rR_FZ: subject correation for full correation global signal regressed
%netmats1R_FZ: subject correation for full correation 
%netmats5R_FZ: subject correation for partial correation

%IN THE FORM: (group, ROI, ROI)
%all_full: group for full correalation with global singal regression
%all_full1: group for full correalation
%all_partial: group for partial correalation

%% start the program 
if corrset==1 %global mean regressed correaltion 
    %get R value to fisher Z trans for mean and sig_corr data 
    netmatsR=nets_netmats(ts,0,'rcorr'); % regress out the global mean --> leads to negative corrleations --> should remove negative corr --> done for mean arival time of blood 
    netmatsZ=nets_netmats(ts,1,'rcorr'); % regress out the global mean --> leads to negative corrleations --> should remove negative corr --> done for mean arival time of blood 
    %THIS SHOULD BE DONE FOR TASK NOT RESTING STATE FMRI ANALYSIS     
    %Macey, Paul M., et al. "A method for removal of global effects from fMRI time series." Neuroimage 22.1 (2004): 360-366.
elseif corrset==2 %partial correaltion 
    netmatsZ=nets_netmats(ts,1,'ridgep');
    netmatsR=nets_netmats(ts,0,'ridgep');
    %THIS SHOULD BE DONE FOR TASK NOT RESTING STATE FMRI ANALYSIS     
    %Macey, Paul M., et al. "A method for removal of global effects from fMRI time series." Neuroimage 22.1 (2004): 360-366.
elseif corrset==3 %for full correaltion 
    netmatsZ=nets_netmats(ts,1,'corr'); %  this is used for resting state,
    netmatsR=nets_netmats(ts,0,'corr'); %  this is used for resting state,
elseif corrset==0
       %FOR EACH SUBJECT
        starttime=1;
       for subji=1:ts.Nsubjects %ALL SUBJECTS USE SAME ROI
            %get ts timepoints for each subject
            endtime=starttime+ts.NtimepointsPerSubject-1;
            %run PPI on only one subject
            onesubj_tsts=ts.ts(starttime:endtime,:);
            starttime=endtime+1;
           [Onsetfile, Movementfile,RT,dt,fMRI_T0,HParam] = SPM_files_PPI(timmingfile,TR,fmri_length);
           save_name='PPIdata';
           [PPI_beta,PPI_tfac,PPI_trRV]=wb_gppi(onesubj_tsts,Onsetfile, Movementfile, save_name,RT,dt,fMRI_T0,HParam);
           %PPI_beta for this will need to be size (ROI*ROI)--in vector
           %form          
           netmatsZ(subji,:)=PPI_beta; %not sure if in same order here? 
            clear PPI_beta
            clear PPI_tfac
            clear PPI_trRV
       end       
end

%fisher Z trasnform these R values 
%netmatsR_FZ=0.5*log((1+netmatsR)./(1-netmatsR)); %why do this out here? 

 %get significnat conections 
[p_uncorrectedZ,p_correctedZ]=nets_glm_grp(netmatsZ,design,contrast,grp,0); 
 
%get mean corr values for group 
lengthgroup=length(group);
for groupi=1:lengthgroup %want for each subject 
    %get the group means for netmats values     
    [ZnetG,MnetG]=nets_groupmean(netmatsZ(group{groupi},:),0);%,lengthtrail);   % test whichever netmat you're interested in; returns Z values from one-group t-test and group-mean netmat
    all_full(groupi,:,:)=MnetG;
    all_full_Z(groupi,:,:)=ZnetG;

    %group_maps=sprintf('%s/cluster_IM_summary_slices',ica_file);
    %% get nets hierarch and nets web --> topolgy mesaurments leave as option? 
  %  if groupi==1
   %     [hier_order,linkages]=nets_hierarchy(Znet1rG,Znet5G,ts,group_maps); 
  %  end
    %netsweb 
   % netwebdir=sprintf('netsweb_group_%d',groupi);
  %  nets_netweb(Znet1rG,Znet5G,ts,group_maps,netwebdir)  
   %saveas(gcf,'noZ_control_ped','jpg')
end

%get the map of mean corr for only significant values 
sig_thres=0.95;%alpha of 0.05 for significnce 
   %for full correaltion 
[sig_corr_all] =sig_corr_task(all_full,p_correctedZ,sig_thres,design,contrast);
[sig_uncorr_all]=sig_corr_task(all_full,p_uncorrectedZ,sig_thres,design,contrast);
    
%from paper:
%http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0032766
%full correlation seems to be the most reliable for global/short coenctions
%

end



