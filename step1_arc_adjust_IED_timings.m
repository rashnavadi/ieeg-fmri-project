%% Written by Tahereh, December 31, 2024

% Define the parent directory for the IED .txt files
parent_dir = '/work/levan_lab/eegfmri_epilepsy/';
output_dir_base = '/work/levan_lab/Tara/'; % Base output directory

% Load the segmentation configuration file
load('/work/levan_lab/Tara/scripts/segments.mat', 'segments');

% Loop through each subject
for subj_num = 1:70
    % Construct the subject ID
    subject_id = sprintf('ICE%03d', subj_num);
    fprintf('Processing Subject: %s\n', subject_id);
    
    % Check if subject has segment data
    if ~isfield(segments, subject_id)
        fprintf('No segment data for subject: %s\n', subject_id);
        continue;
    end
    
    % Define subject-specific directories
    ied_dir = fullfile(parent_dir, subject_id, '3_EEG', '3_Events', 'IED');
    output_dir = fullfile(output_dir_base, ['Subject_', subject_id]);
    
    % Ensure the output directory exists
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end


    % Get all potential IED files
    all_files = dir(fullfile(ied_dir, '*_IED*.txt'));

    % Filter files to include only those with '_IED[0-9]' and exclude 'exclude'
    ied_files = all_files(arrayfun(@(f) ...
        ~isempty(regexp(f.name, '_IED[0-9]+\.txt$', 'once')) && ...
        isempty(regexp(f.name, 'exclude', 'once')), all_files));

    % Loop through the filtered files
    for ied_file_idx = 1:length(ied_files)
        ied_file = ied_files(ied_file_idx).name;
        fprintf('Processing file: %s\n', ied_file);
        ied_path = fullfile(ied_dir, ied_file);

        % Extract the run name (e.g., Run1a, Run2a, etc.) from the file name
        run_name = regexp(ied_file, 'Run\d+[a-z]?', 'match', 'once');
        if isempty(run_name)
            fprintf('Could not extract run name from file: %s\n', ied_file);
            continue;
        end

        % Check if run segmentation data exists
        if ~isfield(segments.(subject_id), run_name)
            fprintf('No segment data for %s (Run: %s)\n', subject_id, run_name);
            continue;
        end

        % Get the start and end volumes for the segment
        segment_range = segments.(subject_id).(run_name);
        segment_start = segment_range(1);
        segment_end = segment_range(2);

        % Calculate the segment's time offset
        segment_offset_time = segment_start * 1.5; % 1.5s TR
        
        fprintf('Processing %s (Run: %s, Start: %d, End: %d)\n', ...
            ied_file, run_name, segment_start, segment_end);
        
        % Read the IED timings (only the first column)
        file_id = fopen(ied_path, 'r'); % Open the file
        file_data = textscan(file_id, '%f%*s%*s', 'Delimiter', '\t', 'TreatAsEmpty', 'NIL');
        fclose(file_id);

        % Extract the numeric IED timings
        ied_timings = file_data{1}; % This will skip "NIL" and treat them as empty

        % Debug: Check IED timings
        disp('IED Timings (First Column Only):');
        disp(ied_timings);
        
        % Adjust timings relative to the segment start
        adjusted_timings = ied_timings + segment_offset_time; % Apply the offset in seconds

        % Save the adjusted timings to a new file
        output_file = fullfile(output_dir, [ied_file(1:end-4), '_adjusted.txt']);
        writematrix(adjusted_timings, output_file, 'Delimiter', ' ');

        fprintf('Adjusted timings saved: %s\n', output_file);
    end
end
