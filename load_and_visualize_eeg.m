%% Written by Tahereh, December 30, 2024

% Script Name: load_and_visualize_eeg.m
% Description: Load EEG .bin file, IED timing .txt file, and visualize a segment of the EEG.

% Define EEG directory
% EEG_DIR = '/work/levan_lab/eegfmri_epilepsy/';
EEG_DIR = '/Users/trashnavadi/Documents/Data_Analysis/2024/postdoc/Levan/analysis/arc_cluster/3_EEG/2_Cleaned/ICE062_Run1_Cleaned/IED_Cleaned';

% Parameters (Update with your file paths)
binFilePath = fullfile(EEG_DIR, 'ICE062_Run1_IED_3000000Samples_55C_2500Hz.bin'); % Path to the .bin file
iedTimingFile = fullfile(EEG_DIR, 'ICE062_Run1a_IED1.txt'); % Path to the .txt file with IED timings
numChannels = 55; % Number of EEG channels
sampleRate = 2500; % Sampling rate in Hz
channelIdx = 1; % Channel to visualize

% Load the .bin file
fileID = fopen(binFilePath, 'r');
eegData = fread(fileID, inf, 'single'); % Assuming float32 binary data
fclose(fileID);

% Reshape the data to [channels x samples]
numSamples = length(eegData) / numChannels; % Number of samples per channel
eegData = reshape(eegData, [numChannels, numSamples]);

% Load IED timings
iedTimings = load(iedTimingFile); % Assuming .txt file contains one timing per line

% Plot the first 5 seconds of the selected channel
durationSeconds = 5; % Duration to plot
startSample = 1; % Starting sample index
endSample = startSample + durationSeconds * sampleRate - 1;

% Ensure the range is within bounds
endSample = min(endSample, size(eegData, 2));

% Extract the EEG segment for the selected channel
time = (startSample:endSample) / sampleRate; % Time vector in seconds
eegSegment = eegData(channelIdx, startSample:endSample);

% Plot the EEG data
figure;
plot(time, eegSegment, 'b');
title(sprintf('EEG Data - Channel %d (First %d seconds)', channelIdx, durationSeconds));
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Overlay IED markers
hold on;
iedInRange = iedTimings(iedTimings >= startSample & iedTimings <= endSample); % Find IEDs in range
iedTimes = (iedInRange - startSample) / sampleRate; % Convert to seconds relative to the plot
for i = 1:length(iedTimes)
    xline(iedTimes(i), 'r', 'LineWidth', 1.5); % Red vertical lines for IEDs
end
hold off;
legend('EEG Signal', 'IED Timings');
