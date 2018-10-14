function [group]=group_defn(matcon)
%This finds the group variable based on the .mat file
%EX: group{1}=[1,2,3]; group{2}=[4,5,6];

%copy file name with next .txt extension
unix(sprintf('cp %s.mat %s.mat.txt',matcon,matcon))

%load design matrix based on text extension
matname=sprintf('%s.mat.txt',matcon);
mattxt = importdata(matname);

%get number of groups and number of subjects from mattxt
num_grps=mattxt.data(1);
num_subj=mattxt.data(2);

%reload matrix to read esay with group and subject numbers 
fileID = fopen(matname);
fname='%s';
newname='%f';
for gi=1:(num_grps-1)
    fname=sprintf('%s %s',fname,newname);
end
C = textscan(fileID,fname);
fclose(fileID);

%CONVERT FIRST TO float and shorten 
s1=C{1};
S1short=s1(6:end);
S = sprintf('%s*', S1short{:});
N = sscanf(S, '%f*');

%make one array with all numbers 
group_array(1,:)=N;
for gi=2:num_grps
    Cone=C{gi};
    group_array(gi,:)=Cone(6:end);
end


% now define group varaible for each value 
for gi=1:num_grps  
    I= find(squeeze(group_array(gi,:)));
    group{gi}=I;  
end

end

