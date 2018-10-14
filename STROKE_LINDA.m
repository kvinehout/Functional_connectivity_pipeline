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
