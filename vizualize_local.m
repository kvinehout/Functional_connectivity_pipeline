function [ total_mask_sum ] = vizualize_local(sig_corr_all,good_compi,data_ica,filename)
%this will vizualize the node WITHIN average connectivity 

cd(data_ica)
total_mask_sum=zeros(91,109,91,size(sig_corr_all,1));
for coni=1:size(sig_corr_all,1) %for length of comparision 
    new_mask_image_sum=zeros(91,109,91);
    for compi=1:size(sig_corr_all,2)%lengthcomp
        %this gets value for all conections 
        new_value_sum=sig_corr_all(coni,compi);  
         %load comp mask 
        loadname=sprintf('mask_comp_whole%d.nii.gz',good_compi(compi));
        nii=load_nii(loadname);
        c = struct2cell(nii);
        nii_mask(:,:,:) = c{5}; 
        clear nii
        %for all componets 
        new_nii_mask_sum=nii_mask.*new_value_sum;
        new_mask_image_sum=new_mask_image_sum+new_nii_mask_sum;
        clear new_value_sum
        clear new_nii_mask
        compi
    end      
    %sum for all componets 
    total_mask_sum(:,:,:,coni)=new_mask_image_sum;
    clear new_value_sum
    coni
end
%save data 
%save new subtraction file 
    new_total_mask_nii_sum_all=make_nii(total_mask_sum);
    filenamenii=sprintf('%s_norm_all_within.nii.gz',filename);
    save_nii(new_total_mask_nii_sum_all,filenamenii)
    fmri=1;
    fmri_length=size(sig_corr_all,1);
    [status] = save_nii_header(filenamenii,fmri,fmri_length);
    %open the new file 
    openfile=sprintf('fslview %s &',filenamenii);
    unix(openfile)
    
end

