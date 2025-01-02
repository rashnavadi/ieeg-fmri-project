%% Written by Tahereh, December 31, 2024

% Define the parent directory for the IED .txt files
ied_dir = '/Users/trashnavadi/Documents/Data_Analysis/2024/postdoc/Levan/analysis/arc_cluster/ICE062/3_EEG/2_Cleaned/ICE062_Run1_Cleaned/IED_Cleaned';

% Define the TR (Repetition Time in seconds)
TR = 1.5;

% Define the start volumes for each run segment
segmentStartVolumes = struct( ...
    'Run1a', 351, ...
    'Run2a', 101, ...
    'Run3a', 176 ...
); % Start of each segment (in fMRI volumes)

% Check if the IED directory exists
if ~isfolder(ied_dir)
    error('IED directory not found: %s', ied_dir);
end

% Get all IED timing files in the directory
ied_files = dir(fullfile(ied_dir, '*_IED*.txt')); % Match all IED .txt files

% Ensure there are files to process
if isempty(ied_files)
    error('No IED .txt files found in directory: %s', ied_dir);
end

% Loop through each IED file and adjust timings
for ied_file_idx = 1:length(ied_files)
    ied_file = ied_files(ied_file_idx).name; % File name
    ied_path = fullfile(ied_dir, ied_file);  % Full path to the file

    % Extract the run name (e.g., Run1a, Run2a, etc.) from the file name
    run_name = regexp(ied_file, 'Run\d+[a-z]?', 'match', 'once');

    if isempty(run_name)
        error('Could not extract run name from file: %s', ied_file);
    end

    % Find the segment start volume for the run
    if ~isfield(segmentStartVolumes, run_name)
        error('No start volume defined for run: %s', run_name);
    end

    % Extract the start and end volumes for the segment
    segment_range = segments.(subject_id).(run_name); % [start, end] volumes for the segment
    segment_start = segment_range(1); % Start volume (e.g., 0-based)
    segment_end = segment_range(2);   % End volume

    % Adjust offset calculation for 0-based indexing
    segment_offset_time = segment_start * TR; % Start volume in seconds (no "-1" needed)

    fprintf('Processing %s (Run: %s, Segment Start Volume: %d, End Volume: %d)\n', ...
        ied_file, run_name, segment_start, segment_end);


    % Read the IED timings (only the first column)
    file_id = fopen(ied_path, 'r'); % Open the file
    file_data = textscan(file_id, '%f', 'Delimiter', '\t', 'TreatAsEmpty', 'NIL');
    fclose(file_id);

    % Extract the numeric IED timings
    ied_timings = file_data{1}; % This will skip "NIL" and treat them as empty
    
    % Adjust timings relative to the segment start
    adjusted_timings = ied_timings + segment_offset_time; % Apply the offset in seconds
    
    % Save the adjusted timings to a new file
    output_file = fullfile(ied_dir, [ied_file(1:end-4), '_adjusted.txt']); % Save with _adjusted suffix
    writematrix(adjusted_timings, output_file, 'Delimiter', ' ');
    
    fprintf('Adjusted timings saved: %s\n', output_file);
end
