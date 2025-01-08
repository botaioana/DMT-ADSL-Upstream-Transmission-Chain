clear variables;
close all; 

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

bitsPerSymbol = numBitsPerSubcarrier * numActiveSubcarriers; % Total bits per DMT symbol
dataRate = bitsPerSymbol * fs; % Data rate in bits per second (bps)

% Initialize full time-domain signal with guard intervals
timeSignalWithGI = [];

%% Step 1: Emitter
% Loop over each DMT symbol
for i = 1:numSymbols
    disp(['Processing Symbol ', num2str(i), ' of ', num2str(numSymbols)]);
    
    % bits to modulated symbols
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
    figure; subplot(211);
    stem(abs(ifftInput), 'filled');
    title(['Magnitude of IFFT Input - Symbol ', num2str(i)]);
    xlabel('Subcarrier Index');
    ylabel('Magnitude');
    grid on;
end
    % Perform IFFT (frequency domain to time domain)
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

% Generate a time vector
numSamples = length(timeSignalWithGI); % Total samples in the transmitted signal
t = (0:numSamples-1) / fe; % Time vector in seconds

% Visualize the full signal
%figure; 
subplot(212);
plot(t, real(timeSignalWithGI));
title('Time Signal with Multiple Symbols and Guard Intervals');
xlabel('Time [s]');
ylabel('Amplitude');
grid on;

% Display total number of samples
disp('Total Samples in Time Signal with Guard Intervals:');
disp(length(timeSignalWithGI));

%% Step 2:  ADSL Channel
% ADSL Channel Simulation: k*sqrt(f) response and noise addition
% Parameters for the channel
SNR_dB = 30; % Desired Signal-to-Noise Ratio
filterOrder = 64; % FIR filter order, default practice
k = 1; % Scaling factor for attenuation

% Design FIR filter with k*sqrt(f) attenuation
f = linspace(0, 1, filterOrder/2 + 1); % Normalized frequency (0 to 1)

% Define gain corresponding to attenuation k*sqrt(f)
gain = 1 ./ (1 + k * sqrt(f)); % Properly scaled gain profile

% Design the FIR filter
h = fir2(filterOrder, f, gain);

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
noise = sqrt(noisePower) * (randn(size(filteredSignal))) ; % Real white noise
noisySignal = filteredSignal + noise;

% Plot the noisy signal
figure;
plot(t, real(timeSignalWithGI)); hold on;
plot(t, real(noisySignal));
title('Noisy Signal After Passing Through ADSL Channel');
xlabel('Time [s]');
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

% Analyzing the performance of the FIR filter
powerBefore = mean(abs(timeSignalWithGI).^2);
powerAfter = mean(abs(filteredSignal).^2);
disp(['Power Before: ', num2str(powerBefore)]);
disp(['Power After: ', num2str(powerAfter)]);

%% Step 3: Receiver

% Placeholder for reconstructed data
reconstructedBits = zeros(size(dataBits)); % Pre-allocate for efficiency

% Define channel response
f = (0:N-1) / N; % Normalized frequency
H = 1 ./ (1 + k * sqrt(f)); % Corrected channel frequency response (low-pass behavior)
% Compute impulse response
impulseResponse = ifft(H); % Time-domain impulse response of the channel

% Plot impulse response
figure;
stem(abs(impulseResponse), 'filled'); % Magnitude of impulse response
title('Channel Impulse Response');
xlabel('Sample Index');
ylabel('Magnitude');
grid on;
hold on;

% Mark the guard interval
xline(guardIntervalLength, 'r--', 'Guard Interval');
legend('Impulse Response', 'Guard Interval');

% Check if impulse response fits within the guard interval
significantSamples = impulseResponse(abs(impulseResponse) > 1e-3); % Threshold to identify significant response
if max(significantSamples) <= guardIntervalLength
    disp('Impulse response fits within the guard interval.');
else
    disp('Impulse response exceeds the guard interval.');
end

for i = 1:numSymbols
    % Extract the current DMT symbol
    startIdx = (i - 1) * (N + guardIntervalLength) + 1;
    endIdx = i * (N + guardIntervalLength);
    receivedSymbol = noisySignal(startIdx:endIdx);

    % Remove guard interval
    receivedSymbol = receivedSymbol(guardIntervalLength+1:end);
    
    % Perform FFT to demodulate
    freqSignal = fft(receivedSymbol, N);
   
    % Extract active subcarriers
    activeSubcarriersData = freqSignal(activeSubcarriers + 1).';
    if i == 1
        scatterplot(activeSubcarriersData(:)); % Flatten to a vector
        xline([-1.26, -0.63, 0,0.63, 1.26], '--r'); % Add vertical decision boundaries for 16-QAM
        yline([-1.26, -0.63, 0,0.63, 1.26], '--r'); % Add horizontal decision boundaries
        xlim([-1.3, 1.3]);% to make it look pretty
        ylim([-1.3, 1.3]);
        title('QAM Constellation before equalization- First Symbol');
    end
    % Equalize the signal (Zero-Forcing Equalization)
    epsilon = 1e-12; % Small value to avoid division by zero
    eqH = conj(H(activeSubcarriers + 1)) ./ (abs(H(activeSubcarriers + 1)).^2 + epsilon); % Equalize
    equalizedSubcarriers = eqH .* activeSubcarriersData;

    % Debug scatter plot for first symbol after equalization
    if i == 1
        scatterplot(equalizedSubcarriers(:)); % Flatten to a vector
        hold on;
        xline([-1.26, -0.63, 0,0.63, 1.26], '--r'); % Add vertical decision boundaries for 16-QAM
        yline([-1.26, -0.63, 0,0.63, 1.26], '--r'); % Add horizontal decision boundaries
        xlim([-1.3, 1.3]);% to make it look pretty
        ylim([-1.3, 1.3]);
        title('Equalized QAM Constellation with Decision Boundaries - First Symbol');
        hold off;
    end

    % Demodulate QAM symbols for active subcarriers
    rQamSymbols = qamdemod(equalizedSubcarriers, M, 'UnitAveragePower', true);
    
    % Convert symbols back to bits
    rQamIndices = de2bi(rQamSymbols, numBitsPerSubcarrier, 'left-msb');
    
    % Reshape the bits to match the expected structure
    rDataBits = reshape(rQamIndices.', [], 1);
    
    % Verify reconstructed bits
    expectedBitsPerSymbol = numBitsPerSubcarrier * numActiveSubcarriers;
    if length(rDataBits) ~= expectedBitsPerSymbol
        error(['Mismatch in bit lengths for Symbol ', num2str(i)]);
    end

    % Map reconstructed bits back to their original position
    bitStartIdx = (i - 1) * expectedBitsPerSymbol + 1;
    bitEndIdx = i * expectedBitsPerSymbol;
    reconstructedBits(bitStartIdx:bitEndIdx) = rDataBits;
end

% Display total reconstructed bits
disp(['Total Reconstructed Bits: ', num2str(length(reconstructedBits))]);

% Calculate the number of erroneous bits
numErrors = nnz(dataBits ~= reconstructedBits);

% Calculate BER
BER = numErrors / length(dataBits);

% Display results
disp(['Number of Erroneous Bits: ', num2str(numErrors)]);
disp(['Total Number of Transmitted Bits: ', num2str(length(dataBits))]);
disp(['Bit Error Rate (BER): ', num2str(BER)]);

%% Step 4: Calculate SNR

signalPower = mean(abs(filteredSignal).^2); % Use the filtered signal before adding noise
noisePower = mean(abs(noise).^2); % Use the noise vector added to the signal
c_SNR = signalPower / noisePower; % Linear SNR
c_SNR_dB = 10 * log10(c_SNR); % SNR in decibels
targetNoisePower = signalPower / (10^(SNR_dB / 10));
disp(['Calculated SNR = ',num2str(c_SNR), ' = ', num2str(c_SNR_dB), ' dB. The initiated value was ', num2str(SNR_dB), ' dB.']);
disp(['Calculated Noise Power (P_noise): ', num2str(targetNoisePower)]);