function [status] = registartion_a2_stroke_control(path_name,reg_file,ana_file,fMRI_file,ana_file_length,linda_path)
%Run registration on the subjects 
%INPUTS: 
%subject_name this is the subject/condtions names 
%path_name this is the path to the subject folders
%reg_file this is the standard space .nii file to regisatar
%ana_file, this is a list of the name of anaomical nii file for each subject
%(within the subject folder)
%fMRI_file this is a list of the name of FMRI nii file for each subject
%(within the subject folder)


%OUTPUTS:
%MAKES REGISTARED FILES called standard_%s.nii.gz within the newly made
%_ICA_AROMA folder within the subject folder


                   
%run bet with referance image
ref=reg_file;
%standard image processing
%MNI orientation
fslreorient2stdname=sprintf('fslreorient2std %s.nii.gz %s_or.nii',ref,ref);
unix(fslreorient2stdname);
betname=sprintf('bet %s_or.nii %s_or_brain.nii -m',ref,ref); 
unix(betname);


cd(path_name)
%note indivudal code made for subjects with data colelcted on mutiple
for  subject=1:length(ana_file)
    anat=ana_file{subject};
    fmri=fMRI_file{subject};
    %cd(subject_name{subject})
    pwd
    %if files already exsit for anatomical regitration skip this part     
    if exist(sprintf('%s_CollapsedWarp.nii.gz',anat),'file') == 2 %so the file exists
        %do nothing, skip to fmRI (IF SAME ANATOMICAL USED FOR MANY FMRI 
    else                    
                   % fl=sprintf('flirt -in %s.nii.gz -ref %s.nii.gz -out %s_fl.nii.gz',fmri,anat,fmri);
                   % unix(fl)
                    %IF contorl run ANTS
                    current_subj=ana_file{subject};
                    if current_subj(1)~='S' %for NOT stroke subjects (controls, hemodyalisis)
                     %This finds the normalization files to referance from
                                    %anatomical image i=numb
                                    %iternations,t=model,r=reg
                                     ANTSn=sprintf('antsRegistrationSyN.sh -d 3 -f %s_or_brain.nii.gz -m %s.nii.gz -o %s_to_MNI.nii.gz -j 1',ref,anat,anat);
                                    unix(ANTSn)                    
                    %IF stroke run LINDA/ANTS
                    elseif current_subj(1)=='S'%for only stroke subjects
                                    %RUN LINDA on anatomical data GET MASK
                                    anatshortBR{1}=anat(1:(end-(ana_file_length+1)));
                                    anat_short=anat(1:(end-(ana_file_length+1)));
                                    [lesion_masks ] = STROKE_LINDA(path_name,anatshortBR,linda_path);
                                    cd(path_name)
                                    %ADD MASK HERE 
                                    ANTSn=sprintf('antsRegistrationSyN.sh -d 3 -f %s_or_brain.nii.gz -m %s.nii.gz -o %s_to_MNI.nii.gz -x %s/linda/finalLindaLesion_binv.nii.gz -j 1',ref,anat,anat,anat_short);
                                    unix(ANTSn)             
                    end
                                     %this transforms anat to the reg image 
                                    warp=sprintf('WarpImageMultiTransform 3 %s.nii.gz %s_warpimtransform.nii.gz -R %s_or_brain.nii.gz %s_to_MNI.nii.gz1Warp.nii.gz %s_to_MNI.nii.gz0GenericAffine.mat', anat,anat,ref,anat,anat);
                                    unix(warp)
                                    %collapse the warp file
                                    collapsewarp=sprintf('antsApplyTransforms -d 3 -o [%s_CollapsedWarp.nii.gz,1] -r %s_or_brain.nii.gz -t %s_to_MNI.nii.gz1Warp.nii.gz -t %s_to_MNI.nii.gz0GenericAffine.mat',anat,ref,anat,anat);
                                    unix(collapsewarp)                                 
    end
    
       %% FMRI 
                                     %get the number of frames for FMRI
                                     hislicen=sprintf('PrintHeader %s.nii.gz | grep Dimens | cut -d '','' -f 4 | cut -d '']'' -f 1',fmri);
                                     [s,hislice]=unix(hislicen);
                                     %get the TR for FMRI data 
                                      trn=sprintf('PrintHeader %s.nii.gz | grep "Voxel Spac" | cut -d '','' -f 4 | cut -d '']'' -f 1',fmri);
                                      [s,tr]=unix(trn);
                                      %convert character to number
                                      hislice=str2num(hislice);
                                      tr=str2num(tr);
                                     %finish registration for FMRI
                                     %replicate the warped image into a 4D warped image for
                                    %FMRI
                                    IMdiffmath=sprintf('ImageMath 3 four_d_deformed_fmri_num%d.nii.gz ReplicateDisplacement %s_CollapsedWarp.nii.gz %d %d 0',subject,anat,hislice,tr);
                                    unix(IMdiffmath)
                                    %replicates the template image into 4D
                                    IMrepmath=sprintf('ImageMath 3 four_d_temp_fmri_num%d.nii.gz ReplicateImage %s_or_brain.nii.gz %d %d 0',subject,ref,hislice,tr);
                                    unix(IMrepmath)
                                    %apply the warped data to the 4D data
                                    applyants=sprintf('antsApplyTransforms -d 4 -o %s_standard.nii.gz -t four_d_deformed_fmri_num%d.nii.gz -r four_d_temp_fmri_num%d.nii.gz -i %s.nii.gz',fmri,subject,subject,fmri);
                                    unix(applyants)
                                    %clear variables
                                    clear tr
                                    clear hislice   
               
           
       %% RUN THIS AFTER REGISTRATION 
       %bias correction after registration  
       unix(sprintf('N4BiasFieldCorrection -d 4 -i %s_standard.nii.gz -o %s_standard_bias.nii.gz',fmri,fmri));
       %denoise
       %need to run fls motion correction with .par file for
       %AROMA 
       mcf2=sprintf('mcflirt -in %s_standard_bias.nii.gz -out %s_standard_mcf.nii.gz -plots',fmri,fmri);
       unix(mcf2)
       %need to run melodic AMROMA for each data set
       AROMA=sprintf('python2 ICA_AROMA.py -in %s_standard_bias.nii.gz -out %s_ICA_AROMA -mc %s_standard_mcf.nii.gz.par -m %s_or_brain_mask.nii.gz',fmri,fmri,fmri,ref);
       unix(AROMA)
       %done with trials      
    cd(path_name)
  
end

  status =1;
end


