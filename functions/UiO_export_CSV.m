% Work in progress and under construction. Idea is to have a generalized export function to save CSV files according to previous calculations, parameters...
% benjamin.thuerer@gmail.com

function UiO_export_CSV(data_struct,subj_name,EEG,logFile)

if nargin < 1
    error('provide at least data_struct. See help UiO_preprocessing')
end

% not needed here? Or generalized from entry of CSV
% if isempty(EEG)
%     if str2double(data_struct.load_data) == 0
%         [EEG,logFile] = UiO_load_data(data_struct,subj_name,'ica_cleaned');   
%     else
%         [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
%     end
% end

% I am afraid global variable is the only possibility to merge across alle
% tests / subjects...
exist export_CSV var
if ans == 0
    global export_CSV
    export_CSV.duration = {};
    export_CSV.duration(1,1) = {"dummy"};
end

% if subj_name already in export_CSV extent testing; else create subj
if strcmp(export_CSV.duration{:,1},subj_name)
    export_CSV =  []; %????
else
    export_CSV.duration(end+1,1) = {subj_name};
end