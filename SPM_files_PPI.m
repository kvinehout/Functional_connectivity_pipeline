function [Onsetfile, Movementfile,RT,dt,fMRI_T0,HParam] = SPM_files_PPI(timmingfile,TR,file_length)
%this function finds the nesisary files in order to run the PPI analysis 

%OUTPUT FILES 
%Onsetfile: the standard onsets .mat file describing the conditions
%                  as used in SPM
                   % USE:http://elden.ua.edu/blog/generating-onset-and-duration-mat-file-for-spm-for-fmri-analysis
%       Movementfile: the standard rp*.txt file containing the movement
%                     regressors
%       Outname: name of the file to save the results (betas) in
% 
%    RT = 1.8; %this is also the TR i thinK? 
%    dt = 0.0545; % timesteps in seconds of microtime resolution%
%	fMRI_T0 = 17; % Corrects for slice-timing offset
%    HParam = 128; % High-pass filter cutoff in seconds   

%RT tthis si the TR? 
%dt this is the timesteps of microtime resolution 
%fMRI_T0 is the correction for slice timing offset
%HParam is the high pass filter cutoff in seconds 




%INPUT FILES
%timmingfile this is the -1,0,1 file one point for each TR this is an
%EXCELL FILE 




%% get the timming varaibles defined 


    RT = TR; %TR
    
    %FIX ALL THESE TO STRANDARD VALUE (MIGHT NEED EDITING LATER)
    dt = 0.0545; % timesteps in seconds of microtime resolution%
	fMRI_T0 = 17; % Corrects for slice-timing offset
    HParam = 128; % High-pass filter cutoff in seconds   






%% get the Movementfile file 

%make a zero movmeent file, flirt was done before here so set to all zeros
move=zeros(file_length,6);
filename='movmentfile.txt';
fid = fopen(filename,'w+');
for ii = 1:size(move,1)
    fprintf(fid,'%20.18f\t',move(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
Movementfile=filename;

%% get the Onsetfile
%load the timming file 
num = xlsread(timmingfile);

%find number of unique numbers in timming file 
C = unique(num);
snum=size(C,1);

%define varaibles 
%This is how the output .m file (for one subject) should look like. There are three conditions (2 task, 1 baseline) and for each condition you have a separate sublist of names, onsets, and durations.
names = cell(1,snum);
onsets = cell(1,snum);
durations = cell(1,snum);

%for each unique number in timming file
for i=1:snum
    %define names varaible
    names{i} = sprintf('activity num %d',i);
    %for this number get onset value from another number        
    %this gives index of all number equal to given value
    k = find(num==C(i));
    dx = [NaN diff(k')];
    removeThis = (dx==1);
    k(removeThis) = [];
    onsets{i} = k'; %KEEP IN SCAN UNITS
    %get the durrations for these onset points 
    %for each onset find until not same number
    for ii=1:length(k)
        %start point
        start=k(ii);
        dur=0;
        while num(start)==C(i) && start < length(num)
            dur=dur+1;
            start=start+1;
            if start >= length(num)
                dur=dur+1;
            end
        end
        durall(ii)=dur;
    end
    durations{i} = durall; %KEEP IN SCAN UNITS
    clear durall
    clear k
    clear dx
end

%save files to onset
Onsetfile='Onsetfile.mat';
save('Onsetfile.mat','names', 'durations', 'onsets')


%% perform SPM analysis on this subject to get varaibles needed for PPI analysis 

%dont want to actually do anything to dataset, just get microtimming and
%stuff 





end

