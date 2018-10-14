function [ status ] = save_nii_header(image,fmri,fmri_length,standard_IM)
%This will go with each save nii file for corrrect header
%image is a sring of the nii file that needs to be saved
%this should be the standard image 
%check if the image is FMRI or MRI image 
standard_IM=' /Users/kaleb/Documents/fsl/data/standard/avg152T1_brain.nii.gz';
if fmri==1
    % fslmerge to fmri_length
    fslmergename=sprintf('fslmerge -t %d_image.nii.gz',fmri_length);
    for i=1:fmri_length
        fslmergeadd = standard_IM;
        fslmergename = strcat(fslmergename,fslmergeadd);
    end
    [status,out]=unix( fslmergename);
    %copy new made header file to the new image 
    filename=sprintf('fslcpgeom %d_image.nii.gz %s',fmri_length,image);
    [status,out]=unix(filename);
elseif fmri==0
    nii_header_string=sprintf('fslcpgeom %s %s',standard_IM,image);
    [status,out]=unix(nii_header_string);
end
end
