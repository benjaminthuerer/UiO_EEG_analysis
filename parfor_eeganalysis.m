%EEGANALYSIS Main analysis script. Calls stuff
%and reads stuff.
% 
% for questions:
% b.e.juel@medisin.uio.no
% sevenius.nilsen@gmail.com
% benjamin.thuerer@kit.edu

function [] = parfor_eeganalysis(sbj,csvname,fldnames1,param)

clear

[csvname, csvpath] = uigetfile('*.csv','choose the csv data in private folder (change file format to "All Files")');
csvname(end-3:end) = [];

tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
addpath('functions\')

% please change the path to your eeglab folder and make sure eeglab runs in
% double precision mode
addpath('Z:\Matlab_Scripts\eeglab14_1_2b\');
eeglab
close

%%Calling the read csv file function
param=readcsv(csvname, csvpath);

%option for parallel pool
distcomp.feature( 'LocalUseMpiexec', false );
pool = parpool();

% subject names
fldnames1 = fieldnames(param);

% parallelized loop for each subject or session in csv
parfor sbj = 1:size(fldnames1,1)
    UiO_par_for_analysis(fldnames1, sbj, param);
end

delete(pool);

end

