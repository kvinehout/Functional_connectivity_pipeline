function [mean_all,SD_all,P,sig_corr_dis,stats] = conection_distance(sig_corr_all,xyz_max)
% THIS FUNTION finds the connection distace between ROI defined connections
%distance defined by COG of ROI defined by cluster program 

%INPUTS: 
%sig_corr_all:this are only signifanct connecionts (comaprision, ROI, ROI)
%xyz_max: this is the xyz for each ROI (3, number of ROI) 

%OUTPUTS: 
%mean_all: mean distance for each contrast
%SD_all: std for each contrast
%P: for contrast and contrast +1 the p value 
%sig_corr_dis: distnace in the form: (coni,node1,node2)
%stats: for contrast and contrast +1 the stats

    
    %% start program 
    
    sig_corr_dis=zeros(size(sig_corr_all,1),size(sig_corr_all,2),size(sig_corr_all,3));  
for coni=1:size(sig_corr_all,1)  
    %get connection paris
    [I,J]=find(squeeze(sig_corr_all(coni,:,:)));
    for pair=1:length(I)
        node1=I(pair);
        node2=J(pair);
        %find xyz of this pair
        max1=xyz_max(:,node1);
        max2=xyz_max(:,node2);
        %get corrdiate distance of xyz pair
        dis(pair)=sqrt((max1(1)-max2(1))^2+(max1(2)-max2(2))^2 +(max1(3)-max2(3))^2);
        sig_corr_dis(coni,node1,node2)=dis(pair);
    end
    if length(I)>=1 %if not all zeros 
        %get average for all coni 
        average_value(coni)=median(dis);
         std_all(coni) = std(dis);
        all_distance{coni,:}=dis;
        coni
        clear dis
        clear I
        clear J
    else
        average_value(coni)=0;
        std_all(coni) = 0;
        all_distance{coni,:}=0;  
    end
end


for coni=1:size(sig_corr_all,1)  
    [MUHAT,SIGMAHAT] = normfit(all_distance{coni});
    mean_all(coni)=MUHAT;
    SD_all(coni)=SIGMAHAT;
end
%}

for i=1:(size(sig_corr_all,1)-1)
[ pi,hi,statsi] = ranksum(all_distance{i},all_distance{i+1});
P(i)=pi;
stats{i}=statsi;
clear statsi
clear pi
end

end

