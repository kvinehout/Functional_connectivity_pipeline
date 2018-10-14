function [ R_all_group,P_all_group,RL_all_group,RU_all_group,sig_all_clinic_R] = clinicalmeasures_a2(path_nameCL,subject_name,all_results,clinical_file,results_name,sig_thres,groupST,group_name)
%run stats on resutls and clinical measures

%INPUTS: 
%path_nameCL: the path to the clinical measurement file
%all_results: (connection,subject) gives value of interest for all subject
    %and all values of interet. 
%clinical_file: this is the excell file with clinical meauerments for eachsubject (IN ORDER OF STATS)
%results_name: Nmae for each result metric EX:results_name{1}='LeftM1-rightM1';
%sig_thres: threshold for comparision P vlaue of 0.05 IN NORMAL P VALUETERMS
%group_name: name of groups EX: group_name{1}='Stroke';group_name{2}='Contorl';
%groupST: group data for data that has clincial measuments 
    %matrix belong to what group
    %EX: group{1}=[1,3,5,6,9]; these subjects are from one group in design IN
    %ORDER OF DESIGN MATRIX
    % group{2}=[2,4,7,8,10]; these subjects are from another group IN
    %ORDER OF DESIGN MATRIX
    %group{n} for n number of groups is experimnet
%subject_name: list of all stroke subject folders included in groupST
    
    
%OUTPUTS: 
%R_all_group: R value of correaltion (group,connection #,clinical measurment #)
%P_all_group: P value for correatlion (group,connection #,clinical measurment #)
%RL_all_group: lower 95% CI (group,connection #,clinical measurment #)
%RU_all_group: upper 95%CI (group,connection #,clinical measurment #)
%sig_all_clinic_R: (group, connection, clinical correlation, Rvalue,Pvalue) for all sig correlaitons



%% load clinical data
%CD INTO AREA WITH CLINICAL MEASURMETNS
cd(path_nameCL)
%load clinical measures
clinical_all=xml2struct(clinical_file);
%for each clinical measure get array
cmC=1;
for cm=2:size(clinical_all.Workbook.Worksheet.Table.Row{3}.Cell,2) %assume first column is subject data
    subC=1;
    %for each subject 
    for sub=1:length(subject_name)%16 %need to define subject automatically
       data=clinical_all.Workbook.Worksheet.Table.Row{sub}.Cell{cm}.Data.Text;
       if sub ==1
           clinical_ST_name{cmC}=data;
       elseif sub>1
            clinical_ST(cmC,subC)=str2num(data);
            subC=subC+1;
       end
    end
    cmC=cmC+1;
end


%% run correaltions 
%for ALL of the results
%divide results input data by group 

for gi=1:length(groupST)
    results_group=all_results(:,groupST{gi});  
    results_ST_all{gi}=results_group;
    clear results_group
end

% RUN FOR EACH GROUP: 
for gi=1:length(groupST)
    %for ped
    results_ST_group=results_ST_all{gi};
    for ri=1:size(results_ST_group,1) %for all results 
        %this if only for ONE clinical measures....want for ALL clincal measure
        for ci=1:size(clinical_ST,1)   %for all clinical meaures
            [R,P,RL,RU] = corrcoef(results_ST_group(ri,:),clinical_ST(ci,:));
            %these are on a diaganal, only want ONE VALUE 
            R_all_group(gi,ri,ci)=R(1,2);
            P_all_group(gi,ri,ci)=P(1,2);
            RL_all_group(gi,ri,ci)=RL(1,2);
            RU_all_group(gi,ri,ci)=RU(1,2);
            clear R
            clear P
            clear RL
            clear RU
            clear results_ST_group
        end
    end
end

%replace any NAN values  with zero
nanarray=isnan(R_all_group);
I=find(nanarray);
[x,y,z]=ind2sub(size(nanarray),I);
R_all_group(x,y,z)=0;
P_all_group(x,y,z)=0;

%get all of the correlations that are less then the lesthold for P vlaue
%test if correlation sig NOT zero 
sig_P_all=P_all_group;
sig_P_all(sig_P_all>sig_thres)=3;
sig_P_all(sig_P_all<=sig_thres)=2;
sig_P_all(sig_P_all==2)=1;
sig_P_all(sig_P_all==3)=0;
%sig_P_all(sig_P_all==nan)=0; %check why nan vlaues 
sig_R_all=sig_P_all.*R_all_group;
I=find(sig_R_all);
[s_taskiR,s_riR,s_ciR]=ind2sub(size(sig_R_all),I); 

%for all significct comparisons
for siR=1:length(I)
    sig_all_clinic_R{siR,1}=group_name{s_taskiR(siR)};
    sig_all_clinic_R{siR,2}=results_name{s_riR(siR)};
    sig_all_clinic_R{siR,3}=clinical_ST_name{s_ciR(siR)};
    sig_all_clinic_R{siR,4}=R_all(s_taskiR(siR),s_riR(siR),s_ciR(siR));
    sig_all_clinic_R{siR,5}= P_all(s_taskiR(siR),s_riR(siR),s_ciR(siR));
end
clear I


%% get figures 

% significnat   sig_all_clinic_R
%for all significnat values 
for i=1:size(sig_all_clinic_R,1)
    size(sig_all_clinic_R,1)
    %find the clinical meausre used
     %for all possible clincal values
     for ii=1:size(clinical_ST_name,2)
           clinused= strfind(clinical_ST_name(ii),sig_all_clinic_R(i,3));
           if clinused{:}>0
                 clinusedall=ii; %shoud be ii
           end 
           clear clinused
     end
     clinicplot=clinical_ST(clinusedall,:);
    clear clinusedall
    %find the results measure used 
    %use all subjects
    clinicplotnew=clinicplot; %this changed if one task not use all stroke subejcts
    % for each sig correlation find what group it belongs to 
    %for all groups
    for gi=1:length(group)
        if strcmp(sig_all_clinic_R(i,1),group_name{gi})==1
            tasknum=gi;
        end  
    end  
    %find what connnection     
      for ii=1:size(results_name,2)
           clinused= strfind(results_name(ii),sig_all_clinic_R(i,2));
           if clinused{:}>0      
                connum=ii;
           end 
           clear clinused
      end
    %get data file 
    taskdata=results_ST_all{tasknum};
    conecplot=taskdata(connum,:);
    clear taskdata
    %need to plot 
    plotdata(:,1)= clinicplotnew;
    plotdata(:,2)=conecplot;
    %X is pair of  conectivity and clinical outcome x in 68 X 2 format 
    %clinical_ST
    %results_ST_all
    clinicalname=sig_all_clinic_R(i,3); 
    clinicalname=clinicalname{:};
    clinicalname(ismember(clinicalname,' ,.:;!()-')) = [];
    t1=sig_all_clinic_R(i,1);
    t1=t1{:};
    t1=t1(3:end);
    t2=sig_all_clinic_R(i,2);
    t2=t2{:};
    connecitonname=sprintf('%s%s',t1,t2); 
    connecitonname(ismember(connecitonname,' ,.:;!()-')) = [];
    clear t1
    clear t2
    %make the figure
    corrplot(plotdata,'varNames',{clinicalname;connecitonname})
    clear clinicalname
    clear connecitonname
    %save the figure 
    save_name=sprintf('correlation number sig 0_1 %d',i);
    saveas(gcf,save_name,'png')
    clear plotdata
end

end

