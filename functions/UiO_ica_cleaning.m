% EEG-data processing for EEG-TMS combined
% Consciousness Study Oslo
% 
% [EEG,logFile] = UiO_ica_cleaning(data_struct,subj_name,EEG,logFile)
% 
% data_struct: structure of the csv-file specified for subject and
%               experiment
% EEG: EEG structure of previous function. If empty [] this function will
%       load the 'after_ica' data (if availeble)
% subj_name: subject name according to csvfile
% logFile: logFile of previous function. If empty [] this function will
%       load the 'after_ica' logFile (if availeble)
%
% This function will clean independent components automized (1) or manually
% by visual inspection (2). Manually by inspection is highly recommendet
% for ICs as the automized method works not very good.
% ICA cleaning works on both, continous and epoched data.
% 
% by questions:
% benjamin.thuerer@kit.edu
%
function [EEG,logFile] = UiO_ica_cleaning(data_struct,subj_name,EEG,logFile)

if nargin < 2
    error('provide at least data_struct and subject name. See help UiO_ica_cleaning')
end


% check if EEG structure is provided. If not, load previous data
if isempty(EEG)
    if str2double(data_struct.load_data) == 0
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,'after_ica');   
    else
        [EEG,logFile] = UiO_load_data(data_struct,subj_name,[],'specific_data');
    end
end

% if eeglab_options not saving properly, icaact might be deleted. This must
% be reconstructed now
if isempty(EEG.icaact)
    if ndims(EEG.data)==3
        for i=1:EEG.trials
            EEG.icaact(:,:,i)=(EEG.icaweights*EEG.icasphere)*squeeze(EEG.data(EEG.icachansind,:,i));
        end
    else
        EEG.icaact = (EEG.icaweights*EEG.icasphere)*EEG.data(EEG.icachansind,:);
    end
end



%% new plugin ICLabel

% check if IC-rejection is automatic by ICLabel (1), automatic by ADJUST (2) or manually (0)
if str2double(data_struct.ica_rejection) == 1
    EEG = iclabel(EEG);
    badICs = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) < 0.5);
    goodICs = setdiff(1:size(EEG.icaact,1),badICs);       
elseif str2double(data_struct.ica_rejection) == 2 
    % time period (ms) for epochs in ICA (just for plotting) and adjust output
    % file
    timePeriod = [-500 500];
    out = 'adjust_output_can_be_deleted';

    %ADJUST needs epoched data. Fake epochs if necessary
    if size(EEG.data,3) == 1
        lag=5; %epoch duration (in seconds)
        disp(['Continuous dataset epoched in ' num2str(lag) ' sec long epochs to compute feature statistics over time']);
        % check whether '5sec' events are present
        si=0; 
        for i=1:length(EEG.event) 
            if(strcmp(EEG.event(1,i).type, '5sec')==1) 
                si=1; 
            end
        end
        if si==0 %add events
            ntrials=floor((EEG.xmax-EEG.xmin)/lag);
            nevents=length(EEG.event);
            for index=1:ntrials
                EEG.event(index+nevents).type='5sec';
                EEG.event(index+nevents).latency=1+(index-1)*lag*EEG.srate; %EEG.srate is the sampling frequency
            end

            EEG=eeg_checkset(EEG,'eventconsistency');
        end

        EEGep = pop_epoch( EEG, {  '5sec'  }, [0 lag], 'newname', [EEG.setname '_ep5'] , 'epochinfo', 'yes');
        %         % removing baseline
        %         EEGep = pop_rmbase( EEGep, [0  0]);
        EEGep = eeg_checkset(EEGep);

        % collects ICA data from EEG
        if isempty(EEGep.icaact)
        warning('EEG.icaact missing! Recomputed from EEG.icaweights, EEG.icasphere and EEG.data');
        % Next instruction: see eeg_checkset
            for i=1:EEGep.trials
                EEGep.icaact(:,:,i) = (EEGep.icaweights*EEGep.icasphere)*squeeze(EEGep.data(EEGep.icachansind,:,i));
            end
        end
    else
        EEGep = EEG;
    end
    
    % run ADJUST und mark bad ICs according to the outcome of ADJUST
    [art] = ADJUST(EEGep,out);
    badICs = art;
    goodICs = setdiff(1:size(EEGep.icaact,1),badICs);
    
elseif str2double(data_struct.ica_rejection) == 0     
    EEG = iclabel(EEG);
    goodICs = [];
    badICs = [];
    
    % for every component: plot the component with ERP, power-spectra,
    % power over trials, and topographic plot    
    for Ci = 1:size(EEG.icaact,1)
        
        h = pop_prop_extended(EEG, 0, Ci, NaN, { 'freqrange' [1 80] }, {'on', 'erp', 'on'}, 0, 'ICLabel');
        suptitle('press space to reject and left mouse click to keep the component')
        
        % seperate good and bad ic's by mouse click / space press
        button = waitforbuttonpress; 
        
        if button==0  
            goodICs = [goodICs, Ci];
        else
            badICs = [badICs, Ci];
        end
        close(h)
    end
    
    % plot ICs again and show the actual marking of the IC. 
    % Press space if you changed your mind
    val = [];
    for Ci = 1:size(EEG.icaact,1)
        if find(goodICs==Ci)
            Stigma = 'good';
        else
            Stigma = 'bad';
        end
        
        h = pop_prop_extended(EEG, 0, Ci, NaN, { 'freqrange' [1 80] }, {'on', 'erp', 'on'}, 0, 'ICLabel');
        suptitle({['marked as ' Stigma] ; ...
                ['if you want to change the mark (regardless of direction g-->b; b-->g) press space'] ; ...
                ['if you want to skip this additional loop press: e']})

        button = waitforbuttonpress;

        % press 'e' to break the loop
        val = double(get(gcf,'CurrentCharacter'));
        if val == 101
            break
        end

        if button ~= 0
            if find(badICs==Ci)
                idx_Ci = badICs == Ci;
                badICs(idx_Ci) = [];
                goodICs(end+1) = Ci;
            else
                badICs(end+1) = Ci;
                idx_Ci = goodICs == Ci;
                goodICs(idx_Ci) = [];
            end
        end
        close(h)
    end 
    goodICs = sort(goodICs);
    badICs = sort(badICs);
else
    error('ica_rejection is not 0, 1, or 2 in the csvfile')
end

% reject bad components on single trial or continous data
decompProj = EEG.icawinv(:, goodICs)*eeg_getdatact(EEG, 'component', goodICs, 'reshape', '2d');

if ndims(decompProj) == 3
    decompProj = reshape(decompProj,size(decompProj,1),EEG.pnts,EEG.trials);
    EEG.data(EEG.icachansind,:,:) = decompProj;
else
    EEG.data(EEG.icachansind,:) = decompProj;
end

% convert data to single and delete unnecissary matrices (save disc space)
EEG.icaact  = [];
EEG.icawinv     = EEG.icawinv(:,goodICs);
EEG.icaweights  = EEG.icaweights(goodICs,:);
EEG.specicaact  = [];
EEG.specdata    = [];
EEG.reject      = [];

% remove EOG channels if exist

if find(cell2mat(strfind({EEG.chanlocs.labels},'EOG')))
    label_cell = strfind({EEG.chanlocs.labels},'EOG');
    for i = 1:length(label_cell)
        if ~isempty(label_cell{i} == 1)
            label_idx(i) = 1;
        else
            label_idx(i) = 0;
        end
    end
    
    label_idx = label_idx == 1;
    if ndims(EEG.data)==3
        EEG.data(label_idx,:,:) = [];
    elseif ndims(EEG.data)==2
        EEG.data(label_idx,:) = [];
    else
        error('not 2 and not 3 dimensions. What else could it be?')
    end
    EEG.chanlocs(label_idx) = [];
end

disp([num2str(length(badICs)) ' ICs rejected']);

% clean up a bit
EEG.nbchan = size(EEG.data,1);


% loc file entry
logFile{end+1} = {'ica_cleaned',['data is cleaned by ICA. In total ' num2str(length(badICs)) ...
    ' were removed.']};

if str2double(data_struct.plot_always)==1
    UiO_plots(data_struct,subj_name,EEG,logFile);
end

disp('data ICA cleaning is done')

end
