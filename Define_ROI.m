function [lengthcomp,IC_comp_cluster] = Define_ROI(ica_file,path_name,min_cluster_size,standard_IM )
%Find the cluster analyis for networks to get number of clusters 
%inputs: 
%standard_IM= this is the file used for registartion earlier' /Users/kaleb/Documents/fsl/data/standard/avg152T1_brain.nii.gz';
%ica_file=this is the folder that has the stadnard FC netwroks   
    %cd('/Users/kaleb/Documents/DATA/standard_ICA_BM70')
%path_name this is the path to the ICA folder
%min_cluster_size this is the mimium voxel size for a cluster
%subject_name thi is a list of all subject/conditions input into meldoic 
    %EX: subject_name={'/pathtto/C01','/pathto/C02');
    
%outputs: 
%length_subject_name:This is just the length of subject/condtion 
%lengthcomp: This is the number of ROI made from the cluster analysis 
%IC_comp_cluster: this is a list of the ICA number assoicated with each new
%ROI 


%%cd into ida file
cd(path_name)
cd(ica_file)
%split the melodic into each comp
unix('fslsplit melodic_IC.nii.gz split_melodic_IC -t');

%% for new analysis of splitting up the componets max values 
%load the componets melodic_IC
melodic_tstatnii=load_nii('melodic_IC.nii.gz');
c = struct2cell(melodic_tstatnii);
melodic_tstat = c{5}; %(x,y,z)
comp_length=length(melodic_tstat(1,1,1,:));
% left and right brain division 
xl=length(melodic_tstat(:,1,1,1));   %182
yl=length(melodic_tstat(1,:,1,1)); %218
zl=length(melodic_tstat(1,1,:,1)); %182
comp_length=length(melodic_tstat(1,1,1,:));
%left and right divide include the extra pixel 
%divde brain along X (1:46,46:91)
midpt=round(xl/2); 
%preallocate memeroy with zero
leftmask_tstat=zeros(xl,yl,zl,comp_length);
rightmask_tstat=zeros(xl,yl,zl,comp_length);
for yi=1:yl
    for zi=1:zl
        for xi=1:xl
            %left mask 
            if xi >= midpt
             leftmask_tstat(xi,yi,zi,:)=melodic_tstat(xi,yi,zi,:);
            end
            %right mask 
            if xi < midpt
             rightmask_tstat(xi,yi,zi,:)=melodic_tstat(xi,yi,zi,:);
            end
        end
    end
end

%make nii file FOR EACH COMP 
for taski=1:comp_length
    out_nii_Lp = make_nii(leftmask_tstat(:,:,:,taski));
    filename=sprintf('leftmask_tstat%d.nii.gz',taski);
    save_nii(out_nii_Lp,filename)
    fmri=0;
    fmri_length=0;
    [status] = save_nii_header(filename,fmri,fmri_length,standard_IM);
    out_nii_Rp = make_nii(rightmask_tstat(:,:,:,taski));
    filename=sprintf('rightmask_tstat%d.nii.gz',taski);
    save_nii(out_nii_Rp,filename)
    fmri=0;
    fmri_length=0;
    [status] = save_nii_header(filename,fmri,fmri_length,standard_IM);  
end
%make nii and open nii for right tstat
    out_nii_Rt = make_nii(rightmask_tstat);
    filename='rightmask_tstat.nii.gz';
    save_nii(out_nii_Rt,filename)
    fmri=1;
    fmri_length=length(rightmask_tstat(1,1,1,:));
    clear rightmask_tstat
    [status] = save_nii_header(filename,fmri,fmri_length,standard_IM); 
    nii=load_nii('rightmask_tstat.nii.gz');
    c = struct2cell(nii);
    rightmask_tstat=c{5};
    
%for left Tstat
    out_nii_Rt = make_nii(leftmask_tstat);
    filename='leftmask_tstat.nii.gz';
    save_nii(out_nii_Rt,filename)
    fmri=1;
    fmri_length=length(leftmask_tstat(1,1,1,:));
    clear leftmask_tstat
    [status] = save_nii_header(filename,fmri,fmri_length,standard_IM); 
    nii=load_nii('leftmask_tstat.nii.gz');
    c = struct2cell(nii);
    leftmask_tstat=c{5};

%% cluster analysis 
%defult here is a Z value of 2.58
%run this for all of the componets 
timepoint_count=1;
%neg_count=1;
for compi=1:comp_length
    %for left AND right 
    for LRi=1:2
        if LRi ==1
            %need to set threshold and min cluster size 
            filenameP=sprintf('cluster --in=leftmask_tstat%d.nii.gz --thresh=2.58 --minextent=%d --oindex=cluster_index_1_comp%d >cluster_comp%d.txt',compi,min_cluster_size,compi,compi);
            unix(filenameP);
            %get the total number of clusters in the image
            cluster1name=sprintf('fslstats cluster_index_1_comp%d.nii.gz -R',compi);
            [status,out]=unix(cluster1name);
            out_str=textscan(out,'%f %f','Delimiter','\t');
            cluster_size=out_str{1,2}; %save cluster size 
        elseif LRi==2
            filenameP=sprintf('cluster --in=rightmask_tstat%d.nii.gz --thresh=2.58 --minextent=%d --oindex=cluster_index_2_comp%d >cluster_comp%d.txt',compi,min_cluster_size,compi,compi);
            unix(filenameP); 
            %get the total number of clusters in the image
            cluster2name=sprintf('fslstats cluster_index_2_comp%d.nii.gz -R',compi);
            [status,out]=unix(cluster2name);
            out_str=textscan(out,'%f %f','Delimiter','\t');
            cluster_size=out_str{1,2};
        end
        
        %for each cluster
        for clus_i=1:cluster_size
                %save variable relating IC comp to each cluster 
                IC_comp_cluster(timepoint_count)=compi;
                filename_cluster=sprintf('fslmaths -dt int cluster_index_%d_comp%d -thr %d -uthr %d -bin cluster_mask%d_%d_%d',LRi,compi,clus_i,clus_i,clus_i,LRi,compi);
                unix(filename_cluster)
                filename_mask_cluster=sprintf('cluster_mask%d_%d_%d.nii.gz',clus_i,LRi,compi);
                cluster_mask=load_nii(filename_mask_cluster);
                filename=sprintf('mask_comp_whole%d.nii.gz',timepoint_count); %save cluster mask 
                %use center of mass data as x,y,z locations 
                textfilename=sprintf('cluster_comp%d.txt',compi);  
                textfile=importdata(textfilename);
                %flip textfile.data
                textfile.data=flipud(textfile.data);
                xw=textfile.data(clus_i,7);
                yw=textfile.data(clus_i,8);
                zw=textfile.data(clus_i,9);  
                max_value_whole(:,timepoint_count)=[xw,yw,zw]; %problem with anatomical matching
                %mask_comp_index=make_nii(mask_comp);
                save_nii(cluster_mask,filename);
                %correct the header information in the mask ADDED THIS 
                fmri=0;
                fmri_length=0;
                [status] = save_nii_header(filename,fmri,fmri_length,standard_IM);
                if LRi ==1 
                    masked=squeeze(leftmask_tstat(:,:,:,compi)).*cluster_mask.img;
                elseif LRi ==2
                    masked=squeeze(rightmask_tstat(:,:,:,compi)).*cluster_mask.img;
                end
            %SAVE the masked image for summary slices analysis 
                masked_all(:,:,:,timepoint_count)= masked;
                %find max of masked image
                %avg_masked=mean(masked(:));%neat to get ride of zero values
                [maxV,id_max] = max(masked(:));
                [minV,id_min]=min(masked(:));
                save('problem_masked.mat','masked')
                if minV <= 0 && maxV <= 0 %ID the negative T value compsTHESE ARE ANTI CORRELATED
                    [xm,ym,zm]=ind2sub([xl,yl,zl],id_min);
                    max_value(:,timepoint_count)=[xm,ym,zm]; %this is negative   
                else %THESE IF FOR POSSITVE CLUSTER AREAS!
                    [xm,ym,zm]=ind2sub([xl,yl,zl],id_max);
                    %save the x,y,z max value for each cluster for each compont
                    max_value(:,timepoint_count)=[xm,ym,zm];
                end
                %save the postive cluster index 
                individual_cluster(timepoint_count,:,:,:)=cluster_mask.img;
                %get the size of each cluster
                temp_size(:,:,:)=cluster_mask.img;
                I = find(temp_size(:));
                total_cluster_size(timepoint_count)=length(I);
                 %make the max image
                mask_comp=zeros(xl,yl,zl);
                mask_comp(xm,ym,zm)=1;
                filename=sprintf('mask_comp_%d.nii.gz',timepoint_count);
                mask_comp_index=make_nii(mask_comp);
                save_nii(mask_comp_index,filename);
                fmri=0;
                fmri_length=0;
                [status] = save_nii_header(filename,fmri,fmri_length,standard_IM);       
                  timepoint_count=timepoint_count+1; 
                
        end
    end
end
%save all the masked cluster data
masked_im=make_nii(masked_all);
save_nii(masked_im,'cluster_IM_summary_slices.nii.gz');
mask_im_size=size(masked_im);
mask_im_size
fmri=1;
fmri_length=length(masked_all(1,1,1,:));
fmri_length
[status] = save_nii_header('cluster_IM_summary_slices.nii.gz',fmri,fmri_length,standard_IM);

%make file for the max_value varaible 
%max_value %x,y,z,cluster,comp
save('max_location_whiten_half_brain.mat','max_value_whole')
save('max_location_single_whiten_half_brain.mat','max_value') %xyz of cluster max 
save('individual_cluster_whiten_half_brain.mat','individual_cluster')
save('total_Cluster_size_whiten_half_brain.mat','total_cluster_size') %size of cluster 
save('IC_comp_cluster_half_brain.mat','IC_comp_cluster')  %relationship of cluster and ICA numbers 
lengthcomp=timepoint_count-1;
lengthcomp



 
%% MAKE GM ONLY ROI FOR MASK_COMP_WHOLE NII IMAGE
%run segmentation on stardard image
standard_IM_base=standard_IM(1:(end-7));
unix(sprintf('/usr/local/fsl/bin/fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -g -B -b -o %s %s',standard_IM_base,standard_IM_base));
%combine segmentation for just GM 
unix(sprintf('fslmaths %s_seg_0.nii.gz -add %s_seg_1.nii.gz %s_seg_GM.nii.gz',standard_IM_base,standard_IM_base,standard_IM_base));            
gmROI=sprintf('%s_seg_GM.nii.gz',standard_IM_base);
for RI=1:lengthcomp
    ROI_file=sprintf('/Users/kaleb/Documents/DATA/standard_ICA_BM70/melodic_2mm/mask_comp_whole%d.nii.gz',RI);
    unix(sprintf('fslmaths %s -mul %s mask_comp_whole_GM_%d.nii.gz',gmROI,ROI_file,RI));  
end
                 

end

