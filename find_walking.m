function [wi, steps, cad] = find_walking(vm, fs, min_amp, T, delta, ...
    alpha, beta, step_freq)
% Identify walking periods and its features from a raw accelerometry signal
% collected using wearable devices.
%
% Detailed method description was published in:
% Straczkiewicz M., Huang E., Onnela J.-P., A “one-size-fits-most” walking
% recognition method for smartphones, smartwatches, and wearable
% accelerometers, npj Digital Medicine, 2023.
%
% Inputs:
% vm ~      vector magnitude of raw acceleration signal
% fs ~      sampling frequency of data collection (in Hz, e.g., 10)
% min_amp ~ amplitude threshold (in g (gravitational units), e.g., 0.2)
% fw ~      step frequency range (in Hz or steps per second, e.g., [1.4,
%           2.3])
% T  ~      minimum walking duration (in seconds, e.g., 3)
% delta ~   maximum difference between consecutive peaks (in multiplication
%           of 0.05Hz, e.g., 2)
% alpha ~   maximum ratio between dominant peak below and within step
%           frequency range (e.g., 0.6)
% beta ~    maximum ratio between dominant peak above and within step
%           frequency range (e.g., 2.5)
%
% Outputs:
% wi ~      walking indication
% steps ~   total number of steps calculated from the input signal
% cad ~     temporal walking cadence (steps per second)
% 
% Example:
% find_walking(vm, 10, 0.2, 3, 2, 0.6, 2.5, [1.4, 2.3])
% 
%
% Script author:
% Marcin Straczkiewicz, PhD
% mstraczkiewicz@hsph.harvard.edu; mstraczkiewicz@gmail.com
%
% Last modification: 20/12/2022

% vectorize the input
vm = vm(:);

% preallocate memory
wi      = zeros(size(vm));
cad     = zeros(numel(vm)/fs,1);
steps   = 0;

% identify valid seconds based on min_amp
valid   = true(size(vm));
pp      = peak2peak(reshape(vm, [fs length(vm)/fs]))';
valid(repelem(pp, fs) < min_amp) = false;

% trim vm to valid segments only
vm = vm(valid);

% decimate valid
valid = valid(1:fs:end);

if numel(vm) < T*fs  % skip processing (valid signal too short)
    return
end

vm_len = length(vm);

% smooth signal at its ends to remove coin of influence
w = tukeywin(numel(vm), 0.02);
vm = vm .* w;
vm = [zeros(5 * fs, 1); vm; zeros(5 * fs, 1)];

if numel(vm) < 50 * fs  
    vm = repmat(vm, [10 1]);% it's a trick done to process even short signals
    [Cima,freqs] = cwt(vm, fs, 'morse', ...
        'WaveletParameters', [3 90], ...
        'VoicesPerOctave', 48, ...
        'NumOctaves', 4);
    Cima = Cima(:, 5 * fs + 1: 5 * fs + vm_len);
else
    [Cima,freqs] = cwt(vm, fs, 'morse', ...
        'WaveletParameters', [3 90], ...
        'VoicesPerOctave', 48, ...
        'NumOctaves', 4);
    Cima = Cima(:, 5 * fs + 1: end - 5 * fs);
end
% get wavelet coefficients
Cabs = abs(Cima).^2; clear Cima

% interpolate wavelet coefficients over new frequency domain
freqs_linspace = round(min(freqs), 1):0.05:round(max(freqs), 1);
Cabs   = interp2(1:vm_len, freqs, Cabs, ...
    1:vm_len, freqs_linspace', ...
    'linear'); % method

% get location of steps frequency boundaries
[~, loc1] = min(abs(freqs_linspace - step_freq(1)));
[~, loc2] = min(abs(freqs_linspace - step_freq(2)));

% find signification peaks in wavelet coefficients and mark them as ones
D = zeros(size(Cabs, 1), size(Cabs, 2) / fs);
j = 0;
for i = 1:size(Cabs, 2)/fs
    x = zeros(size(Cabs, 1), 1);
    % identify 1s of a signal
    j = j + 1;
    vm_1s_start  = (i - 1) * fs + 1;
    vm_1s_finish = (i - 1) * fs + fs;
    
    % find its peaks
    [pks_locs, pks] = peakseek(sum(Cabs(:, vm_1s_start:vm_1s_finish), 2));
    [~, I] = sort(pks, 'descend');
    pks_locs = pks_locs(I);
    pks = pks(I);

    % are there any peaks inside the step frequency range
    step_pks_locs = find(pks_locs >= loc1 & pks_locs <= loc2, 1, 'first');
    if isempty(step_pks_locs)
        continue
    end

    if pks_locs(1) > loc2
        % highest peak is above step frequency range (it might be running
        % or it might be still walking with higher harmonics)
        if pks(1) / pks(step_pks_locs(1)) < beta
            x(pks_locs(step_pks_locs(1))) = 1;
        end
    elseif pks_locs(1) < loc1
        % highest peak is below step frequency range (might be due to hand swinging)
        if pks(1) / pks(step_pks_locs(1)) < alpha
            x(pks_locs(step_pks_locs(1))) = 1;
        end
    else % highest peak is within step frequency range
        x(pks_locs(step_pks_locs(1))) = 1;                   
    end
    % store in matrix
    D(:,j) = x;
end

% align peaks with valid data
E = zeros(size(D, 1), numel(valid));
E(:, valid) = D;

if T == 1
    e = sum(E, 1);
elseif T > 1 % check periodicity of local maxima
    B = find_continuous_dominant_peaks(E, T, delta);
    e = sum(B, 1);
end

% if there is periodicity of a given duration mark it as ones
wi = zeros(size(e));
wi(e > 0) = 1;

% stretch vector to match original signal
wi = repelem(wi, fs);

% estimate temporal cadence based on step frequency of each second
cad = zeros(size(e));
for i = 1:numel(e)
    ind_freqs = find(B(:, i));
    if ~isempty(ind_freqs) && numel(ind_freqs) == 1
        cad(i) = freqs_linspace(ind_freqs);
    end
end
% sum all temporal cadences to get total number of steps
steps = sum(cad);

