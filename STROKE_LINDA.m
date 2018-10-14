function[ lesion_masks ] = STROKE_LINDA(path_name,stroke_ana_file,linda_path)
%This run linda, we will JUST use the lesion mask here NOT resgistation or
%anything else
%linda path is path to LINDA EX:
%linda_path='/Users/kaleb/Desktop/TEST/STROKE_anatomical_AIM1/LINDA_v0.2.7';
cd(sprintf('%spennTemplate/',linda_path))
%load standard mask 
niiS=load_nii('templateBrainMask.nii.gz');
cs = struct2cell(niiS);
MNI_mask = cs{5}; %(x,y,z) %this is one if brain, zero is no brain
%cd(path_name)   WANT TO BE IN LINDA PATH
%for each stroke subject
for subji=1:length(stroke_ana_file)
    %run this part in R outside of matlab 
    %MAKE R SCRIPT
       cd(sprintf('%s%s',path_name,stroke_ana_file{subji}))%cd into subject folder
       fid = fopen('StrokeLINDA.R','wt');
       %define anatomical data
       Name_t1=sprintf('t1="%s%s/ana_or_brain.nii.gz"',path_name,stroke_ana_file{subji});
       %run linda (can not change referance image?)
       ST_LINDA=sprintf('source(''%slinda_predict_noBET.R'')',linda_path);%,stroke_ana_file{subji}); %R should be open at this point
       fprintf(fid, '%s\n%s\n%s',Name_t1,ST_LINDA);
       %close R scipt creation
       fclose(fid);
       %run R script 
       system(sprintf('R CMD BATCH StrokeLINDA.R StrokeLindaOUT.txt'))
        %this will make lesion without the brain outline
        unix(sprintf('fslmaths %s%s/linda/Mask.lesion3.nii.gz -binv %s%s/linda/Mask.lesion3_binv.nii.gz',path_name,stroke_ana_file{subji},path_name,stroke_ana_file{subji}));
        unix(sprintf('fslmaths %s%s/linda/Mask.lesion3_binv.nii.gz -mas %s%s/linda/BrainMask.nii.gz %s%s/linda/finalLindaLesion',path_name,stroke_ana_file{subji},path_name,stroke_ana_file{subji},path_name,stroke_ana_file{subji}));
        unix(sprintf('fslmaths %s%s/linda/finalLindaLesion.nii.gz -binv %s%s/linda/finalLindaLesion_binv',path_name,stroke_ana_file{subji},path_name,stroke_ana_file{subji}));
   %filename=load_nii('Lesion_in_MNI.nii.gz');
   %save the lesion mask to a string
   %lesion_masks{subji}=sprintf('%s/linda/%s',stroke_ana_file{subji},filename); %lesion in MNI space... this is what we want     
end
lesion_masks=1;
end




 %or make R script with just two lines of code and change this script
    %subject name and path name and LINDA path 
    %load the inverted mask 
    %{
    %THIS IS ONLY IF WE WANT LEsion in subject space 
    cd(sprintf('%s/linda',stroke_ana_file{subji}))
    niiL=load_nii('Lesion_in_MNI.nii.gz');
    cd ..
    cd ..
    cl = struct2cell(niiL);
    subj_lesion = cl{5}; %(x,y,z)
    %make zero arrray,%if mask + lesion =0, set to 0,     %if mask =1, lesion =1, set to 0
    subj_lesion_mask=zeros(size(subj_lesion,1),size(subj_lesion,2),size(subj_lesion,3));
    size(subj_lesion)
    size(MNI_mask)
    subj_both=subj_lesion+MNI_mask; %now have 0, adn 2 as bad, and 1 as good 
    I=find(subj_both==1);%find just lesion area
    [x,y,z]=ind2sub(size(subj_both),I); 
    subj_lesion_mask(x,y,z)=1; %set lesion area to one 
    filename=sprintf('%s_MNI_lesion_mask_ones',stroke_ana_file{subji});
    subj_nii_mask=make_nii(subj_lesion_mask);
    save_nii(subj_nii_mask,filename)
    standard_IM=sprintf('%s/linda/Lesion_in_MNI.nii.gz',stroke_ana_file{subji});
    fmri=0;
    fmri_length=0;
   [status] = save_nii_header(filename,fmri,fmri_length,standard_IM);%change header file 
    %}
    
  