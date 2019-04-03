% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_results(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'after_ica' data (if availeble)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the 'ica_cleaned' logFile (if available)
%
% This function will load a matrix and fill this with the results contained
% in EEG.SWA_res (in a future release maybe generalized EEG.res which is 
% independent of SWA). Afterwards it stores everything in this matrix
% again.
% 
% by questions:
% benjamin.thuerer@kit.edu

function [EEG,logFile] = UiO_results(data_struct,subj_name,EEG,logFile)


if nargin < 1
    error('provide at least data_struct. See help UiO_results')
end


% check if EEG structure is provided. If not, load previous data or
% specific dataset
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'SWA_calculated');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% get results_file provided in CSV
results_file = data_struct.results_file;
[filepath, ~, ext] = fileparts(results_file);
if isempty(ext)
    results_file = [results_file '.mat'];
end

% check if results_file and filepath exist (create if not)
if exist(results_file, 'file') == 0
    if exist(filepath, 'dir') == 0
        mkdir(filepath)
    end
    RES = [];
    header = {};
    header.param = EEG.SWA_header;
    header.subj = containers.Map; % use container dictionary to track sessions
else
    load(results_file);
    load([results_file(1:end-4) '_header.mat']);
end

% add individual results (EEG.SWA_res) to overall RES matrix
% RES is 4D: subject, session, channel, parameter
[ch, par] = size(EEG.SWA_res);
if sum(strcmp(keys(header.subj),subj_name{1})) == 0
    RES(end+1, 1, 1:ch, 1:par) = EEG.SWA_res; 
else
    v = values(header.subj);
    RES(end, v{strcmp(keys(header.subj),subj_name{1})}+1, 1:ch, 1:par) = EEG.SWA_res;
end

if sum(strcmp(keys(header.subj),subj_name{1})) == 0
    header.subj(subj_name{1}) = 1;
else
    header.subj(subj_name{1}) = header.subj(subj_name{1}) + 1;
end

% create header to store subj_name
disp('save results and header')
save(results_file,'RES','-v7.3');
save([results_file(1:end-4) '_header.mat'],'header','-v7.3');

end