function [EEG,logFile] = UiO_SWA(data_struct,subj_name,EEG,logFile)

if nargin < 1
    error('provide at least data_struct. See help UiO_preprocessing')
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


% make sure data is in double
EEG.data = double(EEG.data);

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

[number, duration, ptp_amp, numb_pos, numb_neg, power] = UiO_calc_SWA(EEG, epoch_length, marker, SW_length, SW_amplitude);

results = [number, duration, ptp_amp, numb_pos, numb_neg, power];
results_header = ["number", "duration", "ptp_amp", "numb_pos", "numb_neg", "power"];

EEG.SWA_res = results;
EEG.SWA_header = results_header;

% loc file entry
logFile{end+1} = {'SWA_calculated',['SWA is calculated for an epoch of ' num2str(epoch_length) ...
    ' s.']};
end