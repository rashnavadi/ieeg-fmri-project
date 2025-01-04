%% Written by Tahereh, Jan 1, 2025
% Align IED Timings to Peaks Script with Logging
% Updated to handle mismatched Run numbers in file paths.

% Parameters
window_ms = 500; % Window size in ms (±500 ms)
base_eeg_folder = '/work/levan_lab/eegfmri_epilepsy'; % Root folder for EEG data
base_ied_folder = '/work/levan_lab/Tara/'; % Root folder for IED timings
log_file = fullfile(base_ied_folder, 'step2_processing_log.txt'); % Path to log file

% Open log file for writing
fid_log = fopen(log_file, 'w');
fprintf(fid_log, 'Processing Log - IED Alignment\n');
fprintf(fid_log, 'Date: %s\n\n', datestr(now));
fprintf(fid_log, '------------------------------------------------------------\n');

% Subject list
subject_list = arrayfun(@(x) sprintf('ICE%03d', x), 1:70, 'UniformOutput', false);

% Loop through subjects
for subj_idx = 1:length(subject_list)
    subj_id = subject_list{subj_idx};
    fprintf(fid_log, 'Processing subject: %s\n', subj_id);
    disp(['Processing subject: ', subj_id]);

    % Define subject-specific paths
    subj_eeg_folder = fullfile(base_eeg_folder, subj_id, '3_EEG', '2_Cleaned');
    subj_ied_folder = fullfile(base_ied_folder, subj_id);
    output_folder = fullfile(subj_ied_folder, 'aligned_IED_peaks');

    % Create output folder if it does not exist
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % Find all adjusted IED timing files in the subject folder
    ied_files = dir(fullfile(subj_ied_folder, '*_adjusted.txt'));

    if isempty(ied_files)
        fprintf(fid_log, 'No adjusted IED timing files found for subject: %s\n', subj_id);
        disp(['No adjusted IED timing files found for subject: ', subj_id]);
        continue;
    end

    % Loop through each adjusted IED timing file
    for ied_idx = 1:length(ied_files)
        ied_file_name = ied_files(ied_idx).name;
        ied_file_path = fullfile(subj_ied_folder, ied_file_name);
        fprintf(fid_log, 'Processing IED file: %s\n', ied_file_name);
        disp(['Processing IED file: ', ied_file_name]);

        % Extract numeric run number and IED type from the IED file name
        tokens = regexp(ied_file_name, sprintf('%s_Run(\\d+)[a-z]*_IED(\\d+)_adjusted\\.txt', subj_id), 'tokens');
        if isempty(tokens)
            fprintf(fid_log, 'Could not parse IED file name: %s\n', ied_file_name);
            disp(['Could not parse IED file name: ', ied_file_name]);
            continue;
        end
        numeric_run_number = tokens{1}{1}; % Extract the numeric run value (e.g., '1' from 'Run1a')
        ied_type = tokens{1}{2};

        % Find the corresponding EEG folder
        run_folder = sprintf('%s_Run%s_Cleaned', subj_id, numeric_run_number);
        ied_cleaned_folder = fullfile(subj_eeg_folder, run_folder, 'IED_Cleaned');

        if ~exist(ied_cleaned_folder, 'dir')
            fprintf(fid_log, 'IED_Cleaned folder not found: %s\n', ied_cleaned_folder);
            disp(['IED_Cleaned folder not found: ', ied_cleaned_folder]);
            continue;
        end

        bin_files = dir(fullfile(ied_cleaned_folder, '*.bin'));
        if isempty(bin_files)
            fprintf(fid_log, 'No .bin file found for run: %s\n', run_folder);
            disp(['No .bin file found for run: ', run_folder]);
            continue;
        end

        bin_file_name = bin_files(1).name; % Assuming only one .bin file per run
        bin_file_path = fullfile(ied_cleaned_folder, bin_file_name);
        fprintf(fid_log, 'Processing .bin file: %s\n', bin_file_name);
        disp(['Processing .bin file: ', bin_file_name]);

        % Extract information from the .bin file name
        tokens = regexp(bin_file_name, '_(\d+)Samples_(\d+)C_(\d+)Hz\.bin', 'tokens');
        if isempty(tokens)
            fprintf(fid_log, 'Could not parse .bin file name: %s\n', bin_file_name);
            disp(['Could not parse .bin file name: ', bin_file_name]);
            continue;
        end
        num_samples = str2double(tokens{1}{1});
        num_channels = str2double(tokens{1}{2});
        sampling_rate = str2double(tokens{1}{3});

        % Load the EEG data from the .bin file
        try
            fid_bin = fopen(bin_file_path, 'r');
            eeg_data = fread(fid_bin, [num_channels, num_samples], 'single');
            fclose(fid_bin);
        catch ME
            fprintf(fid_log, 'Error reading .bin file: %s\n', bin_file_name);
            fprintf(fid_log, 'Error message: %s\n', ME.message);
            disp(['Error reading .bin file: ', bin_file_name]);
            continue;
        end

        % Load IED timings
        try
            ied_timings = load(ied_file_path);
        catch ME
            fprintf(fid_log, 'Error reading IED timing file: %s\n', ied_file_name);
            fprintf(fid_log, 'Error message: %s\n', ME.message);
            disp(['Error reading IED timing file: ', ied_file_name]);
            continue;
        end

        % Adjust IED timings to align with peaks
        adjusted_timings = zeros(size(ied_timings));
        for ied = 1:length(ied_timings)
            ied_time = ied_timings(ied);

            % Extract indices for the data window around the IED timing
            window_start = max(1, round((ied_time - window_ms / 1000) * sampling_rate));
            window_end = min(num_samples, round((ied_time + window_ms / 1000) * sampling_rate));

            % Extract the data window
            data_window = eeg_data(1, window_start:window_end);

            % Find the peak within the window
            [~, peak_idx] = max(abs(data_window));
            adjusted_timings(ied) = window_start + peak_idx - 1;
        end

        % Save adjusted timings to a new file
        try
            [~, ied_file_base, ~] = fileparts(ied_file_name);
            adjusted_file_name = [ied_file_base, '_Aligned.txt'];
            adjusted_file_path = fullfile(output_folder, adjusted_file_name);
            save(adjusted_file_path, 'adjusted_timings', '-ascii');
            fprintf(fid_log, 'Adjusted timings saved to: %s\n', adjusted_file_path);
            disp(['Adjusted timings saved to: ', adjusted_file_path]);
        catch ME
            fprintf(fid_log, 'Error saving adjusted timings: %s\n', adjusted_file_name);
            fprintf(fid_log, 'Error message: %s\n', ME.message);
            disp(['Error saving adjusted timings: ', adjusted_file_name]);
        end
    end

    fprintf(fid_log, '------------------------------------------------------------\n');
end

fclose(fid_log);
disp('All subjects processed successfully. Check log file for details.');
