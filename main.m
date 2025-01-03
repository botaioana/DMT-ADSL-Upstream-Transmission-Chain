clear variables;

% Modulation Order
M = 16; % QAM modulation order

% Number of Subcarriers
N = 256; % Total subcarriers per FFT (conforming to ADSL standard)

% Sampling Frequency
fe = 2.208e6; % Sampling frequency for ADSL (N * 4.3125 kHz)

% Guard Interval Fraction
guardIntervalFraction = 1/16; % Fraction of the symbol duration

% Effective Symbol Frequency
fs = fe / (1 + guardIntervalFraction); % Effective symbol frequency (includes guard interval)

% Active Subcarriers
activeSubcarriers = 7:31; % Indices of active subcarriers (tones 7 to 31)

% Number of Symbols
numSymbols = 5; % Number of symbols to simulate

% Total number of bits for all symbols
numBitsPerSubcarrier = log2(M); % Bits per QAM symbol
numActiveSubcarriers = numel(activeSubcarriers); % Total number of active subcarriers
totalBits = numSymbols * numBitsPerSubcarrier * numActiveSubcarriers; % Total number of bits
dataBits = randi([0 1], totalBits, 1); % Generate random binary data

% Initialize full time-domain signal with guard intervals
timeSignalWithGI = [];

%totalCyclicPrefix = zeros(numSymbols * round(N * guardIntervalFraction), 1);
% Loop over each symbol
for i = 1:numSymbols
    disp(['Processing Symbol ', num2str(i), ' of ', num2str(numSymbols)]);

    % Extract bits for the current symbol
    startIdx = (i - 1) * numBitsPerSubcarrier * numActiveSubcarriers + 1;
    endIdx = i * numBitsPerSubcarrier * numActiveSubcarriers;
    symbolBits = dataBits(startIdx:endIdx);

    % Reshape bits into symbols
    dataSymbols = reshape(symbolBits, numBitsPerSubcarrier, []).'; % Group bits into symbols
    qamIndices = bi2de(dataSymbols, 'left-msb'); % Convert to integers
    qamSymbols = qammod(qamIndices, M, 'UnitAveragePower', true); % QAM modulation
if i == 1
    % Display QAM constellation for this symbol
    scatterplot(qamSymbols);
    title(['QAM Modulated Symbols - Symbol ', num2str(i)]);
    xlabel('In-phase');
    ylabel('Quadrature');
    grid on;
end
    % Map QAM symbols to IFFT input
    ifftInput = zeros(N, 1); % Initialize IFFT input
    ifftInput(activeSubcarriers + 1) = qamSymbols; % Assign to active subcarriers
    ifftInput(N/2+2:end) = conj(flipud(ifftInput(2:N/2))); % Ensure symmetry for upper subcarriers
if i == 1
    % Verify IFFT Input
    disp('Size of IFFT Input:');
    disp(size(ifftInput)); % Should be [N, 1]
    isSymmetric = all(ifftInput(N/2+2:end) == conj(flipud(ifftInput(2:N/2))));
    disp(['Hermitian Symmetry Valid: ', num2str(isSymmetric)]);

    % Plot magnitude of IFFT input
    figure;
    stem(abs(ifftInput), 'filled');
    title(['Magnitude of IFFT Input - Symbol ', num2str(i)]);
    xlabel('Subcarrier Index');
    ylabel('Magnitude');
    grid on;
end
    % Perform IFFT
    timeSignal = ifft(ifftInput, N);
if i == 1
    % Check if time-domain signal is real
    isRealSignal = all(abs(imag(timeSignal)) < 1e-10);
    disp(['Time Signal is Real: ', num2str(isRealSignal)]);
end
    % Serialize the time-domain signal
    serializedSignal = reshape(timeSignal, [], 1); % Convert parallel form to a single stream

    % Add Guard Interval
    guardIntervalLength = round(N * guardIntervalFraction);
    cyclicPrefix = serializedSignal(end-guardIntervalLength+1:end); % Extract cyclic prefix
    symbolWithGI = [cyclicPrefix; serializedSignal]; % Add cyclic prefix
if i == 1
    % Verify Cyclic Prefix
    isCyclicPrefixCorrect = all(cyclicPrefix == serializedSignal(end-guardIntervalLength+1:end));
    disp(['Cyclic Prefix is Correct: ', num2str(isCyclicPrefixCorrect)]);
end
    % Concatenate to full time-domain signal
    timeSignalWithGI = [timeSignalWithGI; symbolWithGI];
    %guardStart = i * (N + guardIntervalLength) + 1;
    %guardEnd = guardStart + guardIntervalLength - 1;
    %totalCyclicPrefix(guardStart:guardEnd) = cyclicPrefix;
end

% Visualize the full signal
figure;
plot(real(timeSignalWithGI));
title('Time Signal with Multiple Symbols and Guard Intervals');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Display total number of samples
disp('Total Samples in Time Signal with Guard Intervals:');
disp(length(timeSignalWithGI));

%% ADSL Channel
% ADSL Channel Simulation: k*sqrt(f) response and noise addition
% Parameters for the channel
SNR_dB = 30; % Desired Signal-to-Noise Ratio
filterOrder = 64; % FIR filter order

% Design FIR filter with k*sqrt(f) response
f = linspace(0, 1, filterOrder/2 + 1); % Normalized frequency (0 to 1)
magnitudeResponse = sqrt(f); % Magnitude response (k*sqrt(f))
h = fir2(filterOrder, f, magnitudeResponse); % FIR filter design

% Plot the filter's frequency response
figure;
freqz(h, 1, 512);
title('Frequency Response of ADSL Channel (k\sqrt{f})');
xlabel('Normalized Frequency');
ylabel('Magnitude');

% Filter the transmitted signal through the ADSL channel
filteredSignal = conv(timeSignalWithGI, h, 'same');

% Add AWGN to simulate noise
signalPower = mean(abs(filteredSignal).^2); % Signal power
noisePower = signalPower / (10^(SNR_dB / 10)); % Noise power
noise = sqrt(noisePower) * (randn(size(filteredSignal)) + 1j * randn(size(filteredSignal))); % Complex AWGN
noisySignal = filteredSignal + noise;

% Plot the noisy signal
figure;
plot(real(timeSignalWithGI)); hold on;
plot(real(noisySignal));
title('Noisy Signal After Passing Through ADSL Channel');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Signal before', 'Signal after');
grid on; hold off;

% Display total number of samples in the noisy signal
disp('Total Samples in Noisy Signal:');
disp(length(noisySignal));

% Output Parameters for Verification
disp('Channel Simulation Parameters:');
disp(['SNR (dB): ', num2str(SNR_dB)]);
disp(['Filter Order: ', num2str(filterOrder)]);
