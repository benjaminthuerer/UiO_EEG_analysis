% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_SWA(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'after_ica' data (if availeble)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the 'ica_cleaned' logFile (if available)
%
% This function will define markers from CSV file or marker prvided in EEG 
% structure. Then it will call UiO_calc_SWA where slow-wave parameters 
% and low-freq power will be computed. Finally, 'SWA_res' and 'SWA_header' 
% will be added to the EEG structure containing a matrix with slow-wave
% parameters for each channel.
% 
% for questions:
% benjamin.thuerer@kit.edu

function [EEG,logFile] = UiO_SWA(data_struct,subj_name,EEG,logFile)

if nargin < 1
    error('provide at least data_struct. See help UiO_SWA')
end


% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'ica_cleaned');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end
 
% clean ICA structures for further computation (because EOG was removed)
EEG.icawinv = [];
EEG.icasphere = [];
EEG.icaweights = [];
EEG.icachansind = [];

epoch_length = str2double(data_struct.SW_epoch); %epoch length according to CSV
%Find latency for marker either from real marker or manually from CSV

%% be aware that this works in samples and we are resampling!!!
if str2double(data_struct.SW_marker) == 0
    marker = str2double(data_struct.SW_noMarker);
else
    marker = EEG.times(round(EEG.event(find(strcmp({EEG.event.type},"Wake"))).latency));  % get time in ms for sample of marker
end

if isempty(data_struct.SW_duration)
    SW_length = [250,1000];
elseif str2double(data_struct.SW_duration)==0
    SW_length = [250,1000];
else
    delim_idx = strfind(data_struct.SW_duration,'-');
    SW_length(1) = str2double(data_struct.SW_duration(1:delim_idx-1));
    SW_length(2) = str2double(data_struct.SW_duration(delim_idx+1:end));
end

if isempty(data_struct.SW_amplitude)
    SW_amplitude = 75;
elseif str2double(data_struct.SW_amplitude)==0
    SW_amplitude = 75;
else
    SW_amplitude = str2double(data_struct.SW_amplitude);
end

[EEG, number, duration, ptp_amp, numb_pos, numb_neg, power] = UiO_calc_SWA(EEG, epoch_length, marker, SW_length, SW_amplitude);

results = [number, duration, ptp_amp, numb_pos, numb_neg, power];
results_header = ["number", "duration", "ptp_amp", "numb_pos", "numb_neg", "power"];

EEG.SWA_res = results;
EEG.SWA_header = results_header;

% loc file entry
logFile{end+1} = {'SWA_calculated',['SWA is calculated for an epoch of ' num2str(epoch_length) ...
    ' s.']};
end