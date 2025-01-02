%% Written by Tahereh, Jan 1, 2025

% Align IED Timings to Peaks Script with Logging
% Updated to include a log file for tracking processing progress and errors.

% Parameters
window_ms = 500; % Window size in ms (±500 ms)
base_eeg_folder = '/work/levan_lab/eegfmri_epilepsy'; % Root folder for EEG data
base_ied_folder = '/work/levan_lab/Tara/'; % Root folder for IED timings
log_file = fullfile(base_ied_folder, 'processing_log.txt'); % Path to log file


% Subject list
subject_list = arrayfun(@(x) sprintf('ICE%03d', x), 1:70, 'UniformOutput', false);

% Loop through subjects
for subj_idx = 1:length(subject_list)
    subj_id = subject_list{subj_idx};
    disp(['Processing subject: ', subj_id]);

    % Define subject-specific paths
    subj_eeg_folder = fullfile(base_eeg_folder, subj_id, '3_EEG', '2_Cleaned');
    subj_ied_folder = fullfile(base_ied_folder, ['Subject_', subj_id]);
    output_folder = fullfile(subj_ied_folder, 'aligned_IED_peaks');

    % Create output folder if it does not exist
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % Find all adjusted IED timing files in the subject folder
    ied_files = dir(fullfile(subj_ied_folder, '*_adjusted.txt'));

    if isempty(ied_files)
        disp(['No adjusted IED timing files found for subject: ', subj_id]);
        continue;
    end

    % Loop through each adjusted IED timing file
    for ied_idx = 1:length(ied_files)
        ied_file_name = ied_files(ied_idx).name;
        ied_file_path = fullfile(subj_ied_folder, ied_file_name);
        disp(['Processing IED file: ', ied_file_name]);

        % Extract run information from the IED file name
        tokens = regexp(ied_file_name, sprintf('%s_Run(\\d+[a-z]*)_IED(\\d+)_adjusted\\.txt', subj_id), 'tokens');
        if isempty(tokens)
            disp(['Could not parse IED file name: ', ied_file_name]);
            continue;
        end
        run_number = tokens{1}{1};
        ied_type = tokens{1}{2};

        % Find the corresponding EEG .bin file
        run_folder = sprintf('%s_Run%s_Cleaned', subj_id, run_number);
        ied_cleaned_folder = fullfile(subj_eeg_folder, run_folder, 'IED_Cleaned');

        if ~exist(ied_cleaned_folder, 'dir')
            disp(['IED_Cleaned folder not found: ', ied_cleaned_folder]);
            continue;
        end

        bin_files = dir(fullfile(ied_cleaned_folder, '*.bin'));
        if isempty(bin_files)
            disp(['No .bin file found for run: ', run_folder]);
            continue;
        end

        bin_file_name = bin_files(1).name; % Assuming only one .bin file per run
        bin_file_path = fullfile(ied_cleaned_folder, bin_file_name);
        disp(['Processing .bin file: ', bin_file_name]);

        % Extract information from the .bin file name
        tokens = regexp(bin_file_name, '_(\d+)Samples_(\d+)C_(\d+)Hz\.bin', 'tokens');
        if isempty(tokens)
            disp(['Could not parse .bin file name: ', bin_file_name]);
            continue;
        end
        num_samples = str2double(tokens{1}{1});
        num_channels = str2double(tokens{1}{2});
        sampling_rate = str2double(tokens{1}{3});

        % Load the EEG data from the .bin file
        fid = fopen(bin_file_path, 'r');
        eeg_data = fread(fid, [num_channels, num_samples], 'single');
        fclose(fid);

        % Load IED timings
        ied_timings = load(ied_file_path);

        % Adjust IED timings to align with peaks
        adjusted_timings = zeros(size(ied_timings));
        for ied = 1:length(ied_timings)
            ied_time = ied_timings(ied);

            % Extract indices for the data window around the IED timing
            window_start = round((ied_time - window_ms / 1000) * sampling_rate); % Convert to integer
            window_end = round((ied_time + window_ms / 1000) * sampling_rate); % Convert to integer

            % Ensure the indices are within valid bounds
            if window_start < 1
                window_start = 1;
            end
            if window_end > num_samples
                window_end = num_samples;
            end

            % Extract the data window
            data_window = eeg_data(1, window_start:window_end);

            % Find the peak within the window
            [~, peak_idx] = max(abs(data_window));
            adjusted_timings(ied) = window_start + peak_idx - 1;
        end

        % Save adjusted timings to a new file
        [~, ied_file_base, ~] = fileparts(ied_file_name);
        adjusted_file_name = [ied_file_base, '_Aligned.txt'];
        adjusted_file_path = fullfile(output_folder, adjusted_file_name);
        disp(['Saving adjusted timings to: ', adjusted_file_path]);
        save(adjusted_file_path, 'adjusted_timings', '-ascii');
    end
end

disp('All subjects processed successfully.');
