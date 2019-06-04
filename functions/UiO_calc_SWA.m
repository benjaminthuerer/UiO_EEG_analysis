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
function [EEG, number, duration, ptp_amp, numb_pos, numb_neg, power] = UiO_calc_SWA(EEG, epoch_length, marker, SW_length, SW_amplitude)

% find indices of mastoid channels and reref to the average mastoid
ind_M1 = strcmp({EEG.chanlocs.labels},"TP9");
ind_M2 = strcmp({EEG.chanlocs.labels},"TP10");
ref = mean(EEG.data([find(ind_M1), find(ind_M2)],:),1);
EEG.data = bsxfun(@minus, EEG.data, ref);

% remove M1 and M2 from data
EEG.data([find(ind_M1), find(ind_M2)], :) = [];
EEG.chanlocs([find(ind_M1), find(ind_M2)]) = [];
EEG.nbchan = size(EEG.data,1);
EEG.ref = "avgMastoids";

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
%freqz(sosbp, 2^16, EEG.srate) %plot in Hz not rad/sample
EEG.data =  filtfilt(sosbp, gbp, EEG.data); %apply filter design to data

% epoch data
end_mark = marker;
start_mark = end_mark - epoch_length*1000;
[~,idx_start] = min(abs(EEG.times-start_mark));
[~,idx_end] = min(abs(EEG.times-end_mark));
EEG.data = EEG.data(:,idx_start:idx_end-1);
EEG.times = EEG.times(idx_start:idx_end-1);
EEG.pnts = size(EEG.data,2);

% center data
EEG.data = bsxfun(@minus,EEG.data,mean(EEG.data,2));



%%

% find slow waves which a certain min/max length (SW_length) and amplitude (SW_amplitude) (take zero crossings
% as start and end
slow_waves = {};
number = zeros(size(EEG.data,1),1);
duration = zeros(size(EEG.data,1),1);
ptp_amp = zeros(size(EEG.data,1),1);
numb_pos = zeros(size(EEG.data,1),1);
numb_neg = zeros(size(EEG.data,1),1);
num_iterations = 0;
sign_data = diff(sign(EEG.data),[],2);

for chan_i = 1:size(EEG.data,1)
    slow_waves(chan_i).duration = [];
    slow_waves(chan_i).amplitude = [];
    slow_waves(chan_i).pos_peaks = [];
    slow_waves(chan_i).neg_peaks = [];

    
    %% use sign instaed of while loop
     %check whether diff is working on matrix
    idx_zeros = find(sign_data(chan_i,:) ~=0);
  
    i = 1;
    while i <= length(idx_zeros)-2
        %% diff must be calculated on times not on indices (idx * srate) or (EEG.times(idx_zeros))
        diff_zeros = EEG.times(idx_zeros(i+2)) - EEG.times(idx_zeros(i)); %in ms
        ptp = max(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2))) - min(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2)));
        
        % definition of SW: between 250 and 1000 ms long and peak-to-peak > 75 myV
        if diff_zeros >= SW_length(1) && diff_zeros <= SW_length(2) && ptp > SW_amplitude 
            slow_waves(chan_i).(['sw' num2str(i)]) = EEG.data(chan_i,idx_zeros(i):idx_zeros(i+2));
            slow_waves(chan_i).duration(end+1) = length(slow_waves(chan_i).(['sw' num2str(i)]))*(1000/EEG.srate); % SW duration in ms
            slow_waves(chan_i).amplitude(end+1) = ptp; % SW peak-to-peak amplitute
            
            % check whether slow wave starts with negative or positive part
            if sign(mean(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+1)))) == 1
                % check if more than 3 samples available
                if idx_zeros(i+1)-idx_zeros(i) < 3
                    slow_waves(chan_i).pos_peaks(end+1) = 1;
                elseif idx_zeros(i+2)-idx_zeros(i+1) < 3
                    slow_waves(chan_i).neg_peaks(end+1) = 1;
                else
                    slow_waves(chan_i).pos_peaks(end+1) = length(findpeaks(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+1))));
                    slow_waves(chan_i).neg_peaks(end+1) = length(findpeaks(EEG.data(chan_i,idx_zeros(i+1):idx_zeros(i+2))*-1));
                end
            else
                % check if more than 3 samples available
                if idx_zeros(i+1)-idx_zeros(i) < 3
                    slow_waves(chan_i).pos_peaks(end+1) = 1;
                elseif idx_zeros(i+2)-idx_zeros(i+1) < 3
                    slow_waves(chan_i).neg_peaks(end+1) = 1;
                else
                    slow_waves(chan_i).pos_peaks(end+1) = length(findpeaks(EEG.data(chan_i,idx_zeros(i):idx_zeros(i+1))*-1));
                    slow_waves(chan_i).neg_peaks(end+1) = length(findpeaks(EEG.data(chan_i,idx_zeros(i+1):idx_zeros(i+2))));
                end
            end
            
            num_iterations = num_iterations+1;
            i = i+2;
        else
            i = i+1;
        end
    end
    
    % calculate parameters for each channel
    [number(chan_i)] = length(slow_waves(chan_i).duration) / (epoch_length/60); %number of SW per minute (60sec)!
    [duration(chan_i)] = nanmean(slow_waves(chan_i).duration);
    [ptp_amp(chan_i)] = nanmean(slow_waves(chan_i).amplitude);
    [numb_pos(chan_i)] = nanmean(slow_waves(chan_i).pos_peaks);
    [numb_neg(chan_i)] = nanmean(slow_waves(chan_i).neg_peaks);
end

    
% compute normalized power in db (this requires computed baseline)
path = 'Z:\08_Study_Sleep-deprivation\results\EEG_results';
load([path '_mean_baseline.mat'])
EEG.data = EEG.data.^2;
EEG.data = 10*log10(EEG.data./repmat(baseline',1,size(EEG.data,2)));
norm_method = 'db';
[power] = nanmean(EEG.data,2);

%% compute LZ power
epoch = EEG.srate * 8; %THIS IS HARD CODED : 8 sec!!!!!!!!!!!!!!
MLZ = [];
for i = 1:floor(length(EEG.times)/epoch)
    LZdata = EEG.data(:,(i-1)*epoch+1:i*epoch);
    LZdata = LZdata.^2;
    LZdata = LZdata > median(median(LZdata));
    SSsum = sum(LZdata,2);
    [~,index] = sort(SSsum);
    LZdata = LZdata(index,:);
    
    MLZ(i) = UiO_calc_lz_complexity(LZdata(:), 'exhaustive', 1);
    disp(['LZ: ' num2str(MLZ(i))]);
end

EEG.LZ = mean(MLZ);

end

