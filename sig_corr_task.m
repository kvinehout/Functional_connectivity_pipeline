function   [sig_corr_all] = sig_corr_task(all_full,p_correctedf,sig_thres,design,contrast)
%This finds only the significnat correlations and values for ANY dataset
%with ANY comparisons
%load contrast matrix
%copy file name with next .txt extension
unix(sprintf('cp %s %s.txt',contrast,contrast))
conname=sprintf('%s.txt',contrast);
contxt = importdata(conname);
%load design matrix
%copy file name with next .txt extension
unix(sprintf('cp %s %s.txt',design,design))
%load design matrix based on text extension
matname=sprintf('%s.txt',design);
mattxt = importdata(matname);
%get number of groups and number of subjects from mattxt
num_grps=mattxt.data(1);
num_subj=mattxt.data(2);
NumContrasts=contxt.data(2);

%get a matrix of the SUBJECT DATA... get group data

%get a matrix of the CONTRAST DATA
start=NumContrasts+3+1+num_grps+2; %this is the wrong number,not sure why this value changes with design? 

start=start+16; %not sure why need this? 

contrast_v=zeros(NumContrasts,num_grps);
count=1;
for i=start:(start+NumContrasts-1)
     coni=contxt.textdata{i};
     X = str2num(coni);
     %issue here 
     %Unable to perform assignment because the size of the left side is 1-by-8 and the size of the right side is
        %0-by-0.
     contrast_v(count,:)=X;
     clear coni
    count=count+1;
end


%% if global connectivity 
if size(all_full,3)>1
        %get the group means 
        for gi=1:num_grps
            one_mean=squeeze(all_full(gi,:,:));
            all_mean(gi,:,:)=one_mean;
        end

        %for each contorast
        for ci=1:NumContrasts
            %get differance
            %find positive ones
            [Pin]=find(contrast_v(ci,:) >0);
            %find negative ones
            [Nin]=find(contrast_v(ci,:) <0); 
            %if no negative values 
            if isempty(Nin) == 1
                Diff(ci,:,:)=all_mean(Pin,:,:);
            else
                Diff(ci,:,:)=all_mean(Pin,:,:)-all_mean(Nin,:,:);
            end
            %get the sig theshold corr for each condition
            P_Di=reshape(squeeze(p_correctedf(ci,:)),sqrt(size(p_correctedf,2)),sqrt(size(p_correctedf,2)));
            P_Di(P_Di<sig_thres)=0;
            P_Di(P_Di>=sig_thres)=1;
            %save all data
            sig_corr_all(ci,:,:)=squeeze(Diff(ci,:,:)).*P_Di;    
        end

        %make figure of data
        %figure
        for ci=1:NumContrasts
            %subplot(ci,1,1)
            figure
            imagesc(squeeze(sig_corr_all(ci,:,:)))
            colorbar
            title(contxt.textdata{ci})   
        end
end
%% local connectivity 
if size(all_full,3)==1

        %get the group means 
        for gi=1:num_grps
            one_mean=squeeze(all_full(gi,:));
            all_mean(gi,:)=one_mean;
        end
        %for each contorast
        for ci=1:NumContrasts
            %get differance
            %find positive ones
            [Pin]=find(contrast_v(ci,:) >0);
            %find negative ones
            [Nin]=find(contrast_v(ci,:) <0); 
            %if no negative values 
            if isempty(Nin) == 1
                Diff(ci,:)=all_mean(Pin,:);
            else
                Diff(ci,:)=all_mean(Pin,:)-all_mean(Nin,:);
            end
            %get the sig theshold corr for each condition
            P_Di=squeeze(p_correctedf(ci,:));%,sqrt(size(p_correctedf,2)),sqrt(size(p_correctedf,2)));
            P_Di(P_Di<sig_thres)=0;
            P_Di(P_Di>=sig_thres)=1;
            %save all data
            sig_corr_all(ci,:)=squeeze(Diff(ci,:)).*P_Di;    
        end

        %make figure of data
        %figure
        for ci=1:NumContrasts
            %subplot(ci,1,1)
            figure
            imagesc(squeeze(sig_corr_all(ci,:)))
            colorbar
            title(contxt.textdata{ci})   
        end

end
end

