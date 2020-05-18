function [ nodes,xyz_max,sig_corr_all_large] = visualize_global(location_file,size_file,anatomial_file,sig_corr_all,data_ica,save_file)
%this makes the GC glass brian 
%% for ALL significant conections make node and edge file 
%get to where files are saved
cd(data_ica)
if strcmp(save_file,'network')==1 
cd stats
end
%load data
total_cluster_size=importdata(size_file);
max_value=importdata(location_file); 
all_net_name=importdata(anatomial_file); 
%visuzlaization 
lengthcomp=size(sig_corr_all,2);
nodes=cell(6,lengthcomp);
nodes_num=cell(6,lengthcomp);
for compi=1:lengthcomp
    xyz=max_value(:,compi);  
    x=xyz(1);
    y=xyz(2);
    z=xyz(3);
    %XYZ need to be in mm terms this is for CORRECTED SAVE IMAGE ini 2mm
    %sapce
    x=90-2*x;
    y=-126+2*y;
    z=-72+2*z;
    %save max vlaues
    xyz_max(1,compi)=x;
    xyz_max(2,compi)=y;
    xyz_max(3,compi)=z;
    %save the max value to the nodes file
    sX=sprintf('%d',x);
    nodes{1,compi}=sX;
    nodes_num{1,compi}=sX;
    sY=sprintf('%d',y);
    nodes{2,compi}=sY;
    nodes_num{2,compi}=sY;
    sZ=sprintf('%d',z);
    nodes{3,compi}=sZ;
   nodes_num{3,compi}=sZ;
%save the color (4)..set all to 1? 
    nodes{4,compi}='1';
    nodes_num{4,compi}='1';
%save the cluster size (5)
 cluster_size_temp=total_cluster_size(compi);
 sCST=sprintf('%d',cluster_size_temp);
 nodes{5,compi}=sCST;    
nodes_num{5,compi}=sCST; 
%save the cluster label (6)
nodes{6,compi}=all_net_name{compi}; %maybe set this to 'xxx' for testing
numbers=sprintf('%d',compi);
nodes_num{6,compi}=numbers;
compi
end

%save the nodes into a file 
contrast_node='nodes_ex_file_max_melodic_whole.node';
fid = fopen(contrast_node,'w');
for rowi=1:lengthcomp
fmtString = '%s\t%s\t%s\t%s\t%s\t%s\n';
fprintf(fid,fmtString,nodes{:,rowi});
end
fclose(fid);
%make node file into an exicutable 
node_exicutable=sprintf('chmod +x %s',contrast_node);
unix(node_exicutable)

sig_corr_all_large=zeros(size(sig_corr_all));

for contrasti=1:size(sig_corr_all,1)%[1,2,3,4,5,6,23,24,27,28,7,8,17,18]%
        edge_i=squeeze(sig_corr_all(contrasti,:,:));
        if max(edge_i(:))>0
             %threshold the edge file 
            % edge_i(edge_i<thereshold)=0;    
            %make new edge and node file (%save as .edge and .node 
            contrast_edge=sprintf('contrast%d_.edge',contrasti);
            dlmwrite(contrast_edge,edge_i,' ')
            edge_exicutable=sprintf('chmod +x %s',contrast_edge);
            unix(edge_exicutable)
            savename=sprintf('contrast_%d_axial_NO_color_%s.jpg',contrasti,save_file);
            BrainNet_MapCfg('/Users/kaleb/Documents/fsl/data/standard/MNI152_T1_2mm_brain_or.nii.gz','/Users/kaleb/Documents/Matlab/BrainNetViewer_20181219/Data/SurfTemplate/BrainMesh_Ch2withCerebellum.nv',contrast_node,contrast_edge, '/Users/kaleb/Documents/MU_Grad/dissertation/AIM1/nodde_edge_files/format_bold_axial_nocolor.mat',savename);
            close('BrainNet Viewer')

            clear edge_i
        end
    contrasti
end

end
