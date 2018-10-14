function [fin,fmrilengthALL] = preproc_ana_fmri(ana_file,fmri_file,path_name,removeTR,ROBEX_PATH)
%This transfomrs DCAM to nii for anatomical and fMRI
%inputs:
%removeTR: number of TR to remove
%path_name: path name to subject foldres
%fmri_file full path to FMRI folder from path name
%ana_file path from path name to anatomical data


cd(path_name)
    %FOR ANAATOMICAL DATA
    for anai=1:length(ana_file)
        %for each full file path 
        %anatomical data 
        unix(sprintf('"dcm2niix" -b y -z y -o %s -f "ana" "%s"',ana_file{anai},ana_file{anai}));
        %orient the T1 image
        unix(sprintf('fslreorient2std %s/ana.nii.gz %s/ana_or',ana_file{anai},ana_file{anai}));
        %make T1 folder in subject folder and copy .nii file to this folder
        %(one closest to DWI data colleciton)
        unix(sprintf('mkdir %s/T1',ana_file{anai}));
        unix(sprintf('cp %s/ana_or.nii.gz %s/T1/',ana_file{anai},ana_file{anai}));
        %RUN BIAS CORRECTION ON THE FMRI DATA, this is a normalzation method
        %THIS IS AN ANTS PROGRAM
        N4=sprintf('N4BiasFieldCorrection -d 3 -i %s%s/T1/ana_or.nii.gz -o %s%s/T1/ana_or_bias.nii.gz',path_name,ana_file{anai},path_name,ana_file{anai});
        unix(N4)
        %bet..DONT USE BET USE ROBEX 
          cd(ROBEX_PATH)
          unix(sprintf('./ROBEX %s%s/T1/ana_or_bias.nii.gz %s%s/T1/ana_or_brain.nii.gz',path_name,ana_file{anai},path_name,ana_file{anai}));
         cd(path_name)
    end

    %FOR FMRI DATA
    for fmrii=1:length(fmri_file)
        %for fmri need to remove first 4 TR 
        %anatomical data 
        %output name, input file location
        unix(sprintf('"dcm2niix" -b y -z y -o %s -f "fmri" "%s"',fmri_file{fmrii},fmri_file{fmrii}));        
        %orient the T1 image
        unix(sprintf('fslreorient2std %s/fmri.nii.gz %s/fmri_or',fmri_file{fmrii},fmri_file{fmrii}));
        %realgin the fmri data 
        %cd into fmri file
        cd(fmri_file{fmrii})
        %get length of file
        [status,fmrilength]=unix('fslnvols fmri_or.nii.gz');
        %convert character to number
        fmrilength=str2num(fmrilength);
        %split fmri .nii file
        unix('fslsplit fmri_or.nii.gz splitfmri -t')
        %TO REMOVE TR 
        TREND=fmrilength-1;
        if TREND<10
             unix(sprintf('fslmerge -t shortfmri splitfmri000{%d..%d}.nii.gz',removeTR,TREND));
        elseif TREND<100 && TREND>9
             unix(sprintf('fslmerge -t shortfmri splitfmri000{%d..9}.nii.gz splitfmri00{10..%d}.nii.gz',removeTR,TREND));
        elseif TREND<1000 && TREND>99
              unix(sprintf('fslmerge -t shortfmri splitfmri000{%d..9}.nii.gz splitfmri00{10..99}.nii.gz splitfmri0{100..%d}.nii.gz',removeTR,TREND));
        elseif TREND>999 && TREND<999
            unix(sprintf('fslmerge -t shortfmri splitsplitfmri000{%d..9}.nii.gz splitfmri00{10..99}.nii.gz splitfmri0{100..999}.nii.gz splitfmri{1000..%d}.nii.gz',removeTR,TREND));
        end
        fmrilengthALL(fmrii)=fmrilength;
        %cd back into path 
        cd(path_name) 
         %run mcflirt on the functional image this is like linear registration just among fMRI data
         mcf=sprintf('mcflirt -in %s/shortfmri.nii.gz -out %s/shortfmri_or_mcft.nii.gz',fmri_file{fmrii},fmri_file{fmrii});
         unix(mcf)
         %bet function USE BET, ROBEX DOES NOT WORK FOR FMRI
          betn=sprintf('bet %s/shortfmri_or_mcft.nii.gz %s/shortfmri_or_mcft_brain.nii.gz -F',fmri_file{fmrii},fmri_file{fmrii});
          unix(betn)     
          %RUN BIAS CORRECTION ON THE FMRI DATA, this is a normalzation method
          %THIS IS AN ANTS PROGRAM
           N4=sprintf('N4BiasFieldCorrection -d 4 -i %s/shortfmri_or_mcft_brain.nii.gz -o %s/shortfmri_or_mcft_brain_bias.nii.gz',fmri_file{fmrii},fmri_file{fmrii});
           unix(N4)

    end

fin=1;
    
   
end

