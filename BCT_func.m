function [lambda_all,GE_all,LE_all,Clus_all,Spos_all,Sneg_all,vpos_all,vneg_all,cen_all,r_all,cen_all_edge] = BCT_func( sig_corr_all,threshold )
   
% INPUT:
 %sig_corr_all:p values in form: (contrasts, ROI, ROI) 
 %threshold: threshold to ignore non sig results for calcution in 1-P format (EX: 0.95 is p=0.05)


%OUTPUT: 
%lambda_all:
%GE_all,LE_all: global effiecy for all nodes
%Clus_all: cluster coeffiecnt 
%Spos_all: weighted sum only positive
%Sneg_all: weighted sum with negative
%vpos_all:
%vneg_all:
%cen_all:
%r_all:
%cen_all_edge:

    %threhosld thresholded at r = 0.45 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4322636/
    
%% start program
     
     sig_corr_all(sig_corr_all<threshold)=0;
     
    %get measures with BCT software
    %run for all cases in sig corr 
    %this runs for each subject or contrast
    for i=1:size(sig_corr_all,1) %or maybe only for mean cases? 
        %D=squeeze(sig_corr_dis(i,:,:));
        W=squeeze(sig_corr_all(i,:,:));
        W_pos=W;
        W_pos(W_pos<0)=0;
        %BCT fucntions? 
        %can we get network efficientcy
        %D is the distance matrix 
         [lambda,efficiency,ecc,radius,diameter] = charpath(W_pos);
         %lambda is network path length --> shortest path length (L)
         lambda_all(i)=lambda;
         %local efficiecny
        %gloabl efficiecny
         local=0; 
         GE = efficiency_wei(W_pos,local);
         GE_all(i)=GE;
         %local efficiency
          local=2;
         LE = efficiency_wei(W_pos, local);
         LE_all(i,:)=LE;
         %The average shortest path L of a network is the average of all shortest paths between all pairs of vertices. The diameter of a graph is the longest of all shortest paths. Related to the idea of the average shortest path is that of global efficiency, which is the inverse of the average shortest path
        %can we get Clustering coefficient 
        C=clustering_coef_wu(W_pos);
         Clus_all(i,:)=C;
         G=W_pos; %so looking at weights not distacne, distance is the same for all .. no becuase all not significnce  
         %betweenness centraility
         BC=betweenness_wei(G);
         cen_all(i,:)=BC;
         %Edge betweenness centrality
         [EBC,BC]=edge_betweenness_wei(G); 
         cen_all_edge(i,:)=BC;
         %weighted sum 
         [Spos,Sneg,vpos,vneg] = strengths_und_sign(W);
         Spos_all(i,:)=Spos;
         Sneg_all(i,:)=Sneg;
         vpos_all(i)=vpos;
         vneg_all(i)=vneg;
        
         %Local Assortativity measures the extent to which nodes are connected to
        %   nodes of similar strength (vs. higher or lower strength)
        % [loc_assort_pos,loc_assort_neg] = local_assortativity_wu_sign(W); %gives all NAN?? 
         
         %assortativity coefficient measuremetns (should be negative in
         %stroke and postive and contorls) 
         flag=0;%undirected 
         CIJ=W_pos;
         r = assortativity_wei(CIJ,flag);
         r_all(i)=r;
         
         %{
         % comunity structure for overlapping communites 
         type_clustering='complete';%could be 'single' also 
         M=link_communities(W,type_clustering); %is this output the heriarchy clustering? 
         
         M_all(i)=
         
         %get the moduel
         flag=0;%undirected
         Z=module_degree_zscore(W,M,flag);
         
         %}
         
         i
    end
    
    
    
    
    
    
end

