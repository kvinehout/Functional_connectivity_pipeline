function [fMRI_TIME_all] = feat_hemoR_nohemoR(temp_path,fsfname,path_name,subject_name,ref,fMRI_file,TR,TE,delvol,TIMEFILE,hemo_nohemo,subji,name,BLOCKEVENT)
% THIS FILE RUNS FILTERING STEPS FOR REST DATA
%THIS FILE RUNS FILTERING AND HEMOREGRESSION FOR TASK DATA

%this makes and fun fsf files for each subject that do the follwoing: 
%aply 60s HPF 
%apply whittening --> smooths signal --> high freq attentuation 
%this filteres high pass at 0.1 hertz =10 seconds
%this spacial filters with 5mm gauss



%INNPUTS: 
%temp_path='/Users/kaleb/Documents/MU_Grad/dissertation/AIM2/tmp/'; 
%fsfname:this is the template to use to run fsf either rest or task this is the fsf file
%PATH OF original FSF FILE in the form of design_10_5.fsf (NEEDS VARIBLES NAMES IN FILE)
%path_name path to subject folders
%subject_name subject condtion name folders
%fMRI_TIME number of volumes in fmri files
%ref: is the standar space used durring registration
%fMRI_file is the name of the fMRI_file inside ICA_AROMA to run on
%TR is the TR value for fmri data in sec (ex: 2)
%TE is the TE value for the fmri data in sec (EX: 35)
%delvol is the number of volumes to remove at the begning ofthe file(0:length of file)
%timmingfile:  path and .txt for timming file of zero and ones
%hemo_nohemo set to 1 if perfrom hemo, set to 0 to NOT perform hemo 

%OUTPUTS: makes new data in input files based on the filter FSF file



%run for all subjects
for subj=1:length(subject_name)
    subj
    cd(path_name)
    cd(subject_name{subj})
    
   infile=sprintf('%s.nii.gz',fMRI_file{subj});
    %if file not found skip all this 
    if exist(infile, 'file') == 0
        warning('file not found')
        fMRI_TIME_all(subj)=0;
        %continue
    else
    
    %load the fmri data 
    fmrifilt=load_nii(infile);
    c = struct2cell(fmrifilt);
    melodic_tstat = c{5}; %(x,y,z)
    fMRI_TIME=length(melodic_tstat(1,1,1,:));
    lenx=length(melodic_tstat(:,1,1,1));
    leny=length(melodic_tstat(1,:,1,1));
    lenz=length(melodic_tstat(1,1,:,1));
   %calculate for each fmri dataset 
   TOTAL_VOXELS=fMRI_TIME*lenx*leny*lenz;
   cd(temp_path)
   %copy feat files
   unix(sprintf('cp %s tmpDesign.fsf',fsfname));
   copyname=sprintf('cp tmpDesign.fsf design105%d%d%s.fsf',subj,subji,name); %Make a copy for each run
   unix(copyname)
    %define new variables 
    OUTPUTDIRFILE=sprintf('%s/%s/%s_%s',path_name,subject_name{subj},fMRI_file{subj},name);
    INPUTDATANONII=sprintf('%s/%s/%s',path_name,subject_name{subj},fMRI_file{subj});
   %search and trail number globally
      %replace total voxels
   voxrepl=sprintf('sed -ie ''s#TOTALVOXELS#''%d''#g'' design105%d%d%s.fsf',TOTAL_VOXELS,subj,subji,name);%total number
   unix(voxrepl)
      %replace referance
   refrepl=sprintf('sed -ie ''s#PATHTOSTANDARD#''%s''#g'' design105%d%d%s.fsf',ref,subj,subji,name);%standard nii file WITHOUT .nii.gz
   unix(refrepl)
   %repalce TE
   tename=sprintf('sed -ie ''s#TETIME#''%d''#g'' design105%d%d%s.fsf',TE,subj,subji,name);%EX:35
   unix(tename)
   %repalce TR
   trname=sprintf('sed -ie ''s#TRVALUE#''%d''#g'' design105%d%d%s.fsf',TR,subj,subji,name);%EX: 2.000000
   unix(trname)
   %replace volume number
   timerepl=sprintf('sed -ie ''s#fmriLength#''%d''#g'' design105%d%d%s.fsf',fMRI_TIME,subj,subji,name);%s EX: 300
   unix(timerepl)
   %repalce outputdir
   OUTPUTDIRFILEunix=sprintf('sed -ie ''s#OUTPUTDIRFILE#''%s''#g'' design105%d%d%s.fsf',OUTPUTDIRFILE,subj,subji,name);%this is the output file with number folder
   unix(OUTPUTDIRFILEunix)
    %repalce inputdata
   INPUTDATANONIIunix=sprintf('sed -ie ''s#INPUTDATANONII#''%s''#g'' design105%d%d%s.fsf',INPUTDATANONII,subj,subji,name);%path and file name of input data without .nii 
   unix(INPUTDATANONIIunix)
   %delete first TR 
   delvolunix=sprintf('sed -ie ''s#delvol#''%d''#g'' design105%d%d%s.fsf',delvol,subj,subji,name);
   unix(delvolunix)
   %ONLY IF TASK
   if hemo_nohemo ==1
        %set the timming file with zeros and ones (or negative one) for regression 
        TIMEFILEunix=sprintf('sed -ie ''s#TIMEFILE#''%s''#g'' design105%d%d%s.fsf',TIMEFILE,subj,subji,name);% path and .txt for timming file of zero and ones
        unix(TIMEFILEunix)
        %set 1 or 3 collum text file for block or event related desgin
        BLOCKEVENTunix=sprintf('sed -ie ''s#BLOCKEVENT#''%d''#g'' design105%d%d%s.fsf',BLOCKEVENT,subj,subji,name);% path and .txt for timming file of zero and ones
        unix(BLOCKEVENTunix)
   end
   fMRI_TIME_all(subj)=fMRI_TIME;
   %run feat for each run, this needed for feat for some reason, but not
   %other fsl commands
   featname=sprintf('sh -c ". ${FSLDIR}/etc/fslconf/fsl.sh;${FSLDIR}/bin/feat design105%d%d%s.fsf"',subj,subji,name);
   unix(featname)
   %system(featname)
   
   %same this string to a file
  % fid = fopen(sprintf('feat_%d.txt',subj),'wt');
  %  fprintf(fid, featname);
  %  fclose(fid);

   
%  unix(sprintf('chmod u+x feat_1'))

   %}   
   %RENAME THE OUTPUT FILE so same for rest/task
    cd(path_name)   
   %IF NOHEMO THIS IS THE OUPUT FILE 
   if hemo_nohemo ==1
     %IF HEMO THIS IS THE OUTPUT FILE
     unix(sprintf('cp %s/%s_%s.feat/stats/res4d.nii.gz %s/funccondata_%s.nii.gz',subject_name{subj},fMRI_file{subj},name,subject_name{subj},name))
   else
      %copy final file to subject folder and rename file
      unix(sprintf('cp %s/%s_%s.feat/filtered_func_data.nii.gz %s/funccondata_%s.nii.gz',subject_name{subj},fMRI_file{subj},name,subject_name{subj},name))
   end
  subj
%}
 %fMRI_TIME_all=0;
   % end   
end

end

