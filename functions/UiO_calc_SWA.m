% Calculate slow-wave-activity and define parameters
% @ benjamin.thurer@uio.no

% EEG: EEG data set structure according to EEGlab or preprocessing pipeline
% epoch_length: length of epoch for which SWA should be computed (in s)
% marker: time point of stimulus (in ms). Epoch will be set directly before marker

% details:
% 1: re-reference to linked mastoids
% 2: downsample to 128 hz
% 3: bandpass between 0.5 and 4 Hz (bandstop 0.25 and 10 Hz)
% 4: epoch data
% 5: find slow-waves
% 6: calculate parameters

%% what about SWsaturation?
function [number, timing, duration, ptp_amp, numb_pos, numb_neg] = UiO_calc_SWA(EEG, epoch_length, marker)

% find indices of mastoid channels and reref to the average mastoid
ind_M1 = strcmp({EEG.chanlocs.labels},"TP9");
ind_M2 = strcmp({EEG.chanlocs.labels},"TP10");
EEG.data(find(ind_M1),:) = mean(EEG.data([find(ind_M1), find(ind_M2)],:),1);
EEG = pop_reref(EEG, ind_M1);
EEG.ref = "avgMastoids";

% remove M1 and M2 from data
EEG.data(ind_M2,:) = []; 
EEG.chanlocs(ind_M2) = [];
EEG.nbchan = size(EEG.data,1);

%resample
EEG = pop_resample(EEG, 128, 0.8, 0.4);

%filter using chebyshev type II filter: design lowpass filter with less
%than 3dB of ripple in the passband defined from 0.5 (since this is the
%pass band from the high-pass filter from preprocessing) to 4 Hz
%(stop-bands are 0.25 from preprocessing and 10 Hz)

Rp=3; % 3dB of ripple
Rs=40; % 40 db of  attenuation in the stopband
Wp = 4/(EEG.srate/2); %transition badwith start
Ws = 10/(EEG.srate/2); %stop band

[n, Wn]=cheb2ord(Wp,Ws,Rp,Rs);
[z,p,k]=cheby2(n,Rs,Wn); %design cheby type 2 filter
[sosbp,gbp] = zp2sos(z,p,k);  %create second-order representation (looks better)
freqz(sosbp, 2^16, EEG.srate) %plot in Hz not rad/sample
EEG.data =  filtfilt(sosbp, gbp, EEG.data); %apply filter design to data

% epoch data
end_mark = marker;
start_mark = end_mark - epoch_length*1000;
[~,idx_start] = min(abs(EEG.times-start_mark));
[~,idx_end] = min(abs(EEG.times-end_mark));
EEG.data = EEG.data(:,idx_start:idx_end-1);
EEG.times = EEG.times(idx_start:idx_end-1);

% center data
EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data,2));

% find slow waves which are between 0.25 and 1s long (take zero crossings
% as start and end
slow_waves = {};
number = zeros(size(EEG.data,1),1);
duration = zeros(size(EEG.data,1),1);
ptp_amp = zeros(size(EEG.data,1),1);
numb_pos = zeros(size(EEG.data,1),1);
numb_neg = zeros(size(EEG.data,1),1);
num_iterations = 0;

for chan_i = 1:size(EEG.data,1)
    slow_waves(chan_i).duration = [];
    slow_waves(chan_i).amplitude = [];
    slow_waves(chan_i).pos_peaks = [];
    slow_waves(chan_i).neg_peaks = [];
    [~,idx_zeros] = findpeaks(diff(abs(EEG.data(chan_i,:))));
    i = 1;
    while i <= length(idx_zeros)-2
        %% diff must be calculated on times not on indices (idx * srate) or (EEG.times(idx_zeros))
        diff_zeros = EEG.times(idx_zeros(i+2)) - EEG.times(idx_zeros(i)); %check for ms / s / sampling rate
        ptp = max(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2))) - min(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2))); % check for muV
        if diff_zeros >= 250 && diff_zeros <= 1000 && ptp > 75
            slow_waves(chan_i).(['sw' num2str(i)]) = EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2));
            slow_waves(chan_i).duration(end+1) = length(slow_waves(chan_i).(['sw' num2str(i)]))*(1000/EEG.srate); % SW duration in ms
            slow_waves(chan_i).amplitude(end+1) = ptp; % SW peak-to-peak amplitute
            slow_waves(chan_i).pos_peaks(end+1) = length(findpeaks(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2))));
            slow_waves(chan_i).neg_peaks(end+1) = length(findpeaks(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2))*-1));
            num_iterations = num_iterations+1;
            i = i+2;
        else
            i = i+1;
        end
    end
    
    % calculate parameters for each channel
    [number(chan_i)] = num_iterations / (epoch_length/60); %number of SW per minute (60sec)!
    [duration(chan_i)] = nanmean(slow_waves(chan_i).duration);
    [ptp_amp(chan_i)] = nanmean(slow_waves(chan_i).amplitude);
    [numb_pos(chan_i)] = nanmean(slow_waves(chan_i).pos_peaks);
    [numb_neg(chan_i)] = nanmean(slow_waves(chan_i).neg_peaks);
end
 
end

