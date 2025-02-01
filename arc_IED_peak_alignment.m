% Define base directories
base_dir = '/work/levan_lab/eegfmri_epilepsy';
output_base_dir = '/work/levan_lab/Tara'; % Where shifted IEDs are stored

% base_dir = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/original_ICE/4_Data_and_Analysis/';
% output_base_dir = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/original_ICE/4_Data_and_Analysis/Tara';



% Define the subject list
% subject_order = {'ICE001'};
subject_order = {'ICE001', 'ICE002', 'ICE003', 'ICE004', 'ICE005', 'ICE006', ...
                 'ICE007', 'ICE008', 'ICE009', 'ICE010', 'ICE011', 'ICE012', ...
                 'ICE013', 'ICE014', 'ICE015', 'ICE016', 'ICE017', 'ICE018', ...
                 'ICE019', 'ICE020', 'ICE021', 'ICE022', 'ICE023', 'ICE024', ...
                 'ICE025', 'ICE026', 'ICE027', 'ICE028', 'ICE029', 'ICE030', ...
                 'ICE031', 'ICE032', 'ICE033', 'ICE034', 'ICE035', 'ICE036', ...
                 'ICE037', 'ICE038', 'ICE039', 'ICE040', 'ICE041', 'ICE042', ...
                 'ICE043', 'ICE044', 'ICE045', 'ICE046', 'ICE047', 'ICE048', ...
                 'ICE049', 'ICE050', 'ICE051', 'ICE052', 'ICE053', 'ICE054', ...
                 'ICE055', 'ICE056', 'ICE057', 'ICE058', 'ICE059', 'ICE060', ...
                 'ICE061', 'ICE062', 'ICE063', 'ICE064', 'ICE065', 'ICE066', ...
                 'ICE067', 'ICE068', 'ICE069', 'ICE070'};

% Define a mapping for main channels based on subject and IED type
main_channels_map = containers.Map;
main_channels_map('ICE001_IED1') = {'sLpT1', 'sLpT2', 'sLpT3', 'sLpT4'};
main_channels_map('ICE002_IED1') = {'sLmT4', 'sLmT5'};
main_channels_map('ICE002_IED2') = {'sRmT3', 'sRmT4'};
% ICE003 excluded
main_channels_map('ICE004_IED1') = {'sLmT4', 'sLmT5', 'sLiTcv4', 'sLiTcv5', 'sLiTcv6'};
main_channels_map('ICE005_IED1') = {'sLmT5'};
% ICE006 excluded
main_channels_map('ICE007_IED1') = {'gL35'};
% ICE008 excluded
% ICE009 excluded
main_channels_map('ICE010_IED1') = {'sLaP6', 'sLaP7', 'sLaP8'};
main_channels_map('ICE011_IED1') = {'sLmT5', 'sLmT6'};
main_channels_map('ICE011_IED2') = {'sRFPO5', 'sRFPO6', 'sRFPO7'};
main_channels_map('ICE012_IED1') = {'sRmsubT5'};
main_channels_map('ICE012_IED2') = {'sLlT4', 'sLlT5', 'sLlT6'};
main_channels_map('ICE012_IED3') = {'sLmsubT6'};
main_channels_map('ICE013_IED1') = {'dLH5'};
main_channels_map('ICE013_IED2') = {'gLiT5'};
main_channels_map('ICE013_IED3') = {'sLpT1', 'sLpT2', 'sLpT3', 'sLpT4'};
main_channels_map('ICE014_IED1') = {'dRaIN5'};
main_channels_map('ICE014_IED2') = {'dRmIN1', 'dRmIN2', 'dRmIN3'};
main_channels_map('ICE014_IED3') = {'dLpIN2'};
main_channels_map('ICE014_IED4') = {'dLmT1', 'dLmT2'};
main_channels_map('ICE014_IED5') = {'dRmT2', 'dRmT3', 'dRmT4', 'dRmT5'};
% ICE015 excluded
main_channels_map('ICE016_IED1') = {'dLaH1', 'dLmH1'};
main_channels_map('ICE017_IED1') = {'dRpT1', 'dRT2'};
main_channels_map('ICE017_IED2') = {'dRH3', 'dRH4'};
main_channels_map('ICE018_IED1') = {'dRA1', 'dRA2', 'dRH1', 'dRH2', 'dRH3'};
main_channels_map('ICE018_IED2') = {'dLA1', 'dLA2', 'dLA3', 'dLH1', 'dLH2', 'dLH3', 'dLpH1', 'dLpH2', 'dLpH3'};
main_channels_map('ICE019_IED1') = {'sLipmP8', 'sLspmP8', 'sLpcv8'};
main_channels_map('ICE020_IED1') = {'gLs9', 'gLs10', 'gLi18', 'gLi19', 'gLi20'};
main_channels_map('ICE021_IED1') = {'sLiOF3'};
main_channels_map('ICE021_IED2') = {'sLmT5'};
main_channels_map('ICE022_IED1') = {'LA3', 'LH2', 'LAD2'};
main_channels_map('ICE022_IED2') = {'RA3', 'RAH1', 'RAH2', 'RAH3'};
main_channels_map('ICE023_IED1') = {'dRmsT5', 'dRmsT6', 'dRpsT5', 'dRpsT6'};
main_channels_map('ICE024_IED1') = {'dLpIN1', 'dLpIN2', 'dLpIN3'};
% ICE025 excluded
% ICE026 excluded
main_channels_map('ICE027_IED1') = {'dLasT4', 'dLmsT2', 'dLmsT3'};
main_channels_map('ICE027_IED2') = {'dLA7', 'dLA8'};
main_channels_map('ICE028_IED1') = {'dLaR1', 'dLaR2', 'dLpR4', 'dLpR5', 'dLpR6'};
main_channels_map('ICE028_IED2') = {'dLaR1', 'dLaR2', 'dLpR4', 'dLpR5', 'dLpR6', 'gL5', 'gL6', 'gL7'};
main_channels_map('ICE028_IED3') = {'gL5', 'gL6', 'gL7'};
main_channels_map('ICE029_IED1') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4', 'dRpH1', 'dRpH2', 'dRpH3'};
main_channels_map('ICE029_IED2') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLA2', 'dLA3', 'dLA4'};
main_channels_map('ICE030_IED1') = {'dLmOF1', 'dLaH1', 'dLaH7', 'dLaH8', 'dLpH1', 'dLpH7', 'dLpH8'};
main_channels_map('ICE031_IED1') = {'dLSMA6', 'dLSMA7', 'dLSMA8'};
main_channels_map('ICE031_IED2') = {'dLPM5', 'dLPM6', 'dLPM7'};
main_channels_map('ICE032_IED1') = {'sLpTP3', 'sLpTP4', 'sLpTP5', 'sLpTP6'};
main_channels_map('ICE032_IED2') = {'gLiT28', 'gLiT29', 'gLiT30', 'gLiT31'};
main_channels_map('ICE033_IED1') = {'dRasTg1', 'dRasTg2', 'dRasTg3', 'dRasTg4', 'dRaIN2', 'dRaIN3', 'dRaIN4'};
main_channels_map('ICE033_IED2') = {'dRlOF2', 'dRlOF3', 'dRlOF4', 'dRlOF5', 'dRlOF6', 'dRlOF7', 'dRaIN2', 'dRaIN3', 'dRaIN4'};
main_channels_map('ICE034_IED1') = {'dLA1', 'dLA2', 'dLaH1', 'dLaH2', 'dLaH3'};
main_channels_map('ICE034_IED2') = {'dRaH1', 'dRaH2', 'dRA2', 'dRA3'};
main_channels_map('ICE035_IED1') = {'dRmF3', 'dRmF4', 'dRH5', 'dRH6', 'dRaIN1', 'dRaIN2', 'dRaIN3', 'dRaIN4'};
main_channels_map('ICE036_IED1') = {'dRpIN1', 'dRpIN2', 'dRpIN3', 'dRpIN4', 'dRpIN5', 'dRaIN1', 'dRaIN2', 'dRaIN3', 'dRaIN4', 'dRaIN5', 'dRH1', 'dRH2', 'dRH3'};
main_channels_map('ICE036_IED2') = {'dLH1', 'dLH2', 'dLH3'};
main_channels_map('ICE037_IED1') = {'dRA2', 'dRA3', 'dRA4', 'dRaH1', 'dRaH2', 'dRaH3'};
main_channels_map('ICE038_IED1') = {'dLaH3', 'dLaH4', 'dLaH5', 'dLaH6', 'dLpH1', 'dLpH2', 'dLpH3', 'dLpH4'};
main_channels_map('ICE038_IED2') = {'dLmO7', 'dLmO8', 'dLmO9', 'dLmO10'};
main_channels_map('ICE039_IED1') = {'dLaxH3', 'dLaxH4', 'dLaxH5', 'dLaxH6'};
main_channels_map('ICE039_IED2') = {'dLamTg7', 'dLamTg8', 'dLamTg9', 'dLmmTg7', 'dLmmTg8', 'dLmmTg9', 'dLmmTg10', 'dLpmTg5', 'dLpmTg6', 'dLpmTg7', 'dLpmTg8'};
main_channels_map('ICE040_IED1') = {'dRaH3', 'dRaH4', 'dRaH5', 'dRaH6', 'dRpH1', 'dRpH2', 'dRpH3', 'dRpH4', 'dRpH5'};
main_channels_map('ICE040_IED2') = {'dRpIN2'};
main_channels_map('ICE040_IED3') = {'dLaH3', 'dLA5'};
main_channels_map('ICE041_IED1') = {'dRTcx4', 'dRTcx5', 'dRTcx6', 'dRTP4', 'dRTP5', 'dRTP6'};
main_channels_map('ICE041_IED2') = {'dRmP2', 'dRmP3', 'dRmP4', 'dRmP5', 'dRiP2', 'dRiP3', 'dRiP4', 'dRiP5', 'dRiP6'};
main_channels_map('ICE041_IED3') = {'dRPOP2', 'dRPOP3', 'dRPOP4', 'dRPOP5'};
main_channels_map('ICE041_IED4') = {'dRlOF3', 'dRlOF4', 'dRlOF5', 'dRlOF6'};
main_channels_map('ICE042_IED1') = {'dRpC6', 'dRpC7', 'dRpC8'};
main_channels_map('ICE042_IED2') = {'dRmC1', 'dRmC2', 'dRmC3', 'dRSMA1', 'dRSMA2', 'dRSMA3'};
main_channels_map('ICE042_IED3') = {'dRmC7', 'dRmC8', 'dRaC5', 'dRaC6'};
main_channels_map('ICE042_IED4') = {'dRlOF8', 'dRlOF9', 'dRlOF10'};
main_channels_map('ICE043_IED1') = {'dRaH1', 'dRaH2'};
main_channels_map('ICE043_IED2') = {'dLaH1', 'dLaH2', 'dLA1', 'dLA2', 'dLA3', 'dLpH1', 'dLpH2'};
main_channels_map('ICE044_IED1') = {'dLaH1', 'dLaH2'};
main_channels_map('ICE044_IED2') = {'dRaH1', 'dRaH2', 'dRaH3'};
main_channels_map('ICE044_IED3') = {'dLpIN3', 'dLpIN4'};
main_channels_map('ICE045_IED1') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLpH1', 'dLpH2'};
main_channels_map('ICE045_IED2') = {'dRA1', 'dRA2', 'dRA3', 'dRaH1', 'dRaH2', 'dRaH3'};
main_channels_map('ICE046_IED1') = {'dRUmTg2'};
main_channels_map('ICE046_IED2') = {'dRaIN2'};
main_channels_map('ICE047_IED1') = {'dRaH1', 'dRaH2', 'dRaH3'};
main_channels_map('ICE047_IED2') = {'dLpH1', 'dLpH2', 'dLpH3'};
main_channels_map('ICE048_IED1') = {'dLaH1', 'dLaH2', 'dLpH1', 'dLH2', 'dLpH3'};
main_channels_map('ICE049_IED1') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLA1', 'dLA2', 'dLpH1', 'dLpH2', 'dLpH3'};
main_channels_map('ICE049_IED2') = {'dLA6', 'dLA7', 'dLaH5', 'dLaH6', 'dLaH7', 'dLasT1', 'dLasT2', 'dLasT3', 'dLasT4', 'dLasT5'};
main_channels_map('ICE049_IED3') = {'dLA6', 'dLA7'};
main_channels_map('ICE050_IED1') = {'dLM12', 'dLM13', 'dLM14', 'dLSMA3', 'dLcOP3', 'dLpIN1', 'dLpIN2', 'dLpIN3', 'dLpIN4', 'dLpIN5', 'dLpIN6'};
main_channels_map('ICE051_IED1') = {'dRlmF4', 'dRlmF5', 'dRlmF6', 'dRlmF7', 'dRlmF8', 'dRaIN2', 'dRaIN3'};
main_channels_map('ICE052_IED1') = {'dRpsT1', 'dRpsT2', 'dRpsT3', 'dRpsT4'};
main_channels_map('ICE052_IED2') = {'dRPOP1', 'dRPOP2', 'dRPOP3', 'dRPOP4'};
main_channels_map('ICE053_IED1') = {'dLaH1', 'dLaH2', 'dLaH3', 'dLpH2', 'dLpH3'};
main_channels_map('ICE054_IED1') = {'dRiP5', 'dRiP6'};
main_channels_map('ICE054_IED2') = {'dRpsT5', 'dRpsT6'};
main_channels_map('ICE054_IED3') = {'dRiP5', 'dRiP6', 'dRpsT5', 'dRpsT6'};
main_channels_map('ICE055_IED1') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4'};
main_channels_map('ICE055_IED2') = {'dLaH1', 'dLaH2'};
main_channels_map('ICE056_IED1') = {'dLA3', 'dLaH1', 'dLaH2', 'dLpH1', 'dLpH2'};
main_channels_map('ICE056_IED2') = {'dRA1', 'dRA2', 'dRaH1', 'dRaH2', 'dRpH1', 'dRpH2'};
main_channels_map('ICE057_IED1') = {'dRaIN1', 'dRaIN2', 'dRaIN3', 'dRaIN4'};
main_channels_map('ICE057_IED2') = {'dLmOF5', 'dLmOF6', 'dLmOF7', 'dLmOF8', 'dLmOF9'};
main_channels_map('ICE058_IED1') = {'dLaH2', 'dLaH3', 'dLpH1', 'dLpH2', 'dLA1', 'dLA2', 'dLA3'};
main_channels_map('ICE059_IED1') = {'dRA5', 'dRA6', 'dRA7', 'dRA8', 'dRA9'};
main_channels_map('ICE059_IED2') = {'dLA6', 'dLA7', 'dLA8'};
main_channels_map('ICE059_IED3') = {'dLA2', 'dLA3'};
main_channels_map('ICE059_IED4') = {'dRaH1', 'dRaH2', 'dRA1', 'dRA2'};
main_channels_map('ICE060_IED1') = {'dRA1', 'dRA2', 'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4', 'dRpH1', 'dRpH2', 'dRpH3', 'dRpH4'};
% ICE061 excluded
main_channels_map('ICE062_IED1') = {'dRaH1', 'dRaH2','dRA1', 'dRA2', 'dRA3', 'dRpIN4'};
main_channels_map('ICE062_IED2') = {'dLaH2', 'dLA1', 'dLA2', 'dLpIN1', 'dLpIN2'};
main_channels_map('ICE063_IED1') = {'dLpH1', 'dLpH2', 'dLaH2', 'dLaH3'}; 
main_channels_map('ICE063_IED2') = {'dLA1', 'dLA2'};
main_channels_map('ICE064_IED1') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRaH4'};
main_channels_map('ICE064_IED2') = {'dRsO5', 'dRsO6', 'dRsO7', 'dRsO8'};
main_channels_map('ICE065_IED1') = {'dLpH1', 'dLpH2', 'dLpH3', 'dLaH1', 'dLaH2', 'dLaH3'}; 
main_channels_map('ICE066_IED1') = {'dRpIN1', 'dRpIN2', 'dRpIN3'};
main_channels_map('ICE066_IED2') = {'dRaH1', 'dRaH2', 'dRpH2'};
% ICE067 excluded
% ICE068 excluded
main_channels_map('ICE069_IED1') = {'dLaH1', 'dLaH2', 'dLA1', 'dLA2'};
main_channels_map('ICE069_IED2') = {'dRaH1', 'dRaH2', 'dRaH3', 'dRA1', 'dRA2'};
main_channels_map('ICE070_IED1') = {'dRaIN1', 'dRaIN2'};
main_channels_map('ICE070_IED2') = {'dRasT3', 'dRasT4', 'dRasT5'};
main_channels_map('ICE070_IED3') = {'dRPOP4', 'dRPOP5'};
% Add more subject-IED mappings as needed


% Loop through all subjects
for subj_idx = 1:length(subject_order)
    subject = subject_order{subj_idx};
    fprintf('Processing subject: %s\n', subject);

    output_dir = fullfile(output_base_dir, subject, 'IED_analysis');
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    % Define the new IED directory (where shifted IEDs are stored)
    events_dir = fullfile(output_base_dir, subject);
    
    % List all files in the new events directory before filtering
    all_ied_files = dir(fullfile(events_dir, '*.txt'));
    if isempty(all_ied_files)
        fprintf('No IED files found in directory: %s\n', events_dir);
    else
        fprintf('IED files found:\n');
        for f = 1:length(all_ied_files)
            fprintf('  - %s\n', all_ied_files(f).name);
        end
    end

    % Get all shifted IED files (matching _adjusted.txt format)
    ied_timing_files = dir(fullfile(events_dir, '*_IED*_adjusted.txt')); %% IEDs that have already been shifted/adjusted wrt EEG timings (this IEDs are the original inputs for the rest of the analyses)
    
    % Check if any shifted IED files exist
    if isempty(ied_timing_files)
        fprintf('Skipping subject: %s (No adjusted IED timing files found)\n', subject);
        continue;
    end

    % Get all run folders inside EEG cleaned directory
    eeg_cleaned_dir = fullfile(base_dir, subject, '3_EEG', '2_Cleaned');
    all_entries = dir(fullfile(eeg_cleaned_dir, '*_Cleaned')); % Find all *_Cleaned folders
    run_folders = all_entries([all_entries.isdir]); % Keep only directories
    
    % Process each run folder
    for run_idx = 1:length(run_folders)
        run_folder_name = run_folders(run_idx).name;
     
        % Define full path to EEG cleaned data (new requirement)
        eeg_dir = fullfile(base_dir, subject, '3_EEG', '2_Cleaned', run_folder_name, 'IED_Cleaned');
        
        % Identify EEG file
        eeg_file = dir(fullfile(eeg_dir, '*.bin'));
        if isempty(eeg_file)
            fprintf('Skipping run folder: %s (No EEG file found in %s)\n', run_folder_name, eeg_dir);
            continue;
        end
        eeg_file_path = fullfile(eeg_dir, eeg_file(1).name);
%         fprintf('Processing EEG file: %s\n', eeg_file_path);

        % Process each IED timing file
        for ied_file_idx = 1:length(ied_timing_files)
            ied_timing_file = fullfile(events_dir, ied_timing_files(ied_file_idx).name);
            fprintf('Processing shifted IED timing file: %s\n', ied_timing_file);
            
            % Extract run name from the filename (e.g., ICE049_Run2a_IED3_shifted.txt)
            run_name = regexp(ied_timing_files(ied_file_idx).name, 'Run\d+[a-z]?', 'match', 'once');
            if isempty(run_name)
                fprintf('Could not extract run name from file: %s\n', ied_timing_files(ied_file_idx).name);
                continue;
            end

            % Extract IED type from filename (e.g., ICE002_Run2_IED1.txt → IED1)
            [~, ied_filename, ~] = fileparts(ied_timing_files(ied_file_idx).name);
            ied_filename = strrep(ied_filename, '_adjusted', '');
            ied_parts = split(ied_filename, '_');
            ied_type = ied_parts{end}; % Last part should be IED type (e.g., "IED1")
            
            % Ensure this run matches the cleaned EEG directory
            if ~contains(run_folder_name, run_name)
                fprintf('Skipping %s, as it does not match the expected run: %s\n', ...
                        run_folder_name, run_name);
                continue;
            end
            
            % Print confirmation that we're processing this file
            fprintf('Processing %s (Run: %s)\n', ied_timing_files(ied_file_idx).name, run_name);
                    
            log_file = dir(fullfile(eeg_dir, '*ProcessingLog_&_ChannelLabels.txt'));
            log_file = log_file(~startsWith({log_file.name}, '._')); % Exclude hidden files

            if isempty(log_file)
                continue;
            end
            log_file_path = fullfile(eeg_dir, log_file(1).name);
            fid = fopen(log_file_path, 'r');
            log_contents = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            fclose(fid);
            log_contents = log_contents{1}; % Extract as a cell array of strings
            
            % Locate "Output Parameters:" line
            output_idx = find(contains(log_contents, 'Output Parameters:'));

            % Initialize variables
            sampling_rate = []; n_channels = []; channel_labels = {};

            % Read variables starting from the "Output Parameters" section
            collecting_labels = false;  % Flag to track channel label collection
            for i = (output_idx +1):length(log_contents)
                log_line = strtrim(log_contents{i});

                % Extract variables
                if startsWith(log_line, 's_Rate:')
                    sampling_rate = sscanf(log_line, 's_Rate: %d');
                elseif startsWith(log_line, 'n_Channels:')
                    n_channels = sscanf(log_line, 'n_Channels: %d');
                elseif startsWith(log_line, 'c_Labels (Row):')
                    collecting_labels = true;
                    channel_labels = {}; % Reset labels collection
                elseif collecting_labels
                    % Collect channel labels
                    if isempty(log_line) || contains(log_line, ':')
                        collecting_labels = false;  % Stop collecting if a new section starts
                    else
                        channel_labels = [channel_labels, strsplit(log_line)];
                    end
                end
            end
            % Display extracted information
            fprintf('Sampling Rate: %d Hz\n', sampling_rate);
            fprintf('Number of Channels: %d\n', n_channels);
            fprintf('Channel Labels: %s\n', strjoin(channel_labels, ', '));

            % Validate extraction
            if isempty(sampling_rate) || isempty(n_channels) || isempty(channel_labels)
                return;
            end

           
            % ANALYSES FOR EACH RUN OF EACH SUBJECT FOR DIFFERENT TYPES OF IEDS

            % Define main channel names based on subject and IED type
            key = sprintf('%s_%s', subject, ied_type);
            if isKey(main_channels_map, key)
                main_channel_names = main_channels_map(key);
            else
                main_channel_names = {'DefaultCh1', 'DefaultCh2', 'DefaultCh3'}; % Default
            end
            
            % Read EEG data
            fileID = fopen(eeg_file_path, 'r');
            eeg_data = fread(fileID, [n_channels, Inf], 'float32');
            fclose(fileID);
            n_samples = size(eeg_data, 2);
            
            % Load IED timings
            ied_timings_seconds = load(ied_timing_file);
            ied_timings_samples = round(ied_timings_seconds * sampling_rate);
            
            % Step 1: Find IED peak times using sum of squares and window
            % size of 20 msec
            window_samples = round(0.01 * sampling_rate); % size of +/-10 msec (or 25 samples) was applied to work for subjects with dense IEDs
            % Find main channel indices
            main_channel_indices = find(ismember(channel_labels, main_channel_names));
            disp('Main Channel Names:'); disp(main_channel_names);

            % Initialize shifted timings
            shifted_ied_timings = zeros(size(ied_timings_samples));
            
            % Process each IED
            % Step 1: 
            % calculates the sum of squares of signal values within the ±10 ms window across the main channels
            for i = 1:length(ied_timings_samples)
                ied_sample = ied_timings_samples(i);
                search_start = max(1, ied_sample - window_samples);
                search_end = min(n_samples, ied_sample + window_samples);
                
                % Extract signals and compute sum of squares
                main_signal = eeg_data(main_channel_indices, search_start:search_end);
                main_signal_demeaned = main_signal - mean(main_signal, 2);
                sum_of_squares = sum(main_signal_demeaned.^2, 1);
                
                % Find peak and adjust timing
                [~, max_idx] = max(sum_of_squares);
                shifted_ied_timings(i) = search_start + max_idx - 1;
            end
            
            shifted_ied_timings_seconds = shifted_ied_timings / sampling_rate;
            % Save shifted timings
            writematrix(shifted_ied_timings_seconds, strrep(ied_timing_file, '.txt', '_shifted_peak_timings.txt'), 'Delimiter', 'tab');
            
            % Step 2: Compute average IED waveform (±100 ms)
            averaging_window = round(0.1 * sampling_rate);
            average_ieds = zeros(n_channels, 2 * averaging_window + 1);
            for ch = 1:n_channels
                channel_waveforms = zeros(length(shifted_ied_timings), 2 * averaging_window + 1);
                for i = 1:length(shifted_ied_timings)
                    avg_window_start = max(1, shifted_ied_timings(i) - averaging_window);
                    avg_window_end = min(n_samples, shifted_ied_timings(i) + averaging_window);
                    channel_waveforms(i, :) = eeg_data(ch, avg_window_start:avg_window_end);
                end
                average_ieds(ch, :) = mean(channel_waveforms, 1, 'omitnan');
            end
            writematrix(average_ieds, strrep(ied_timing_file, '.txt', '_average_IEDs.txt'), 'Delimiter', 'tab');
            
            % Step 3: Find peak amplitude within ±50 ms
            small_window = round(0.05 * sampling_rate);
            zero_idx = averaging_window + 1;
            search_start = max(1, zero_idx - small_window);
            search_end = min(size(average_ieds, 2), zero_idx + small_window);
            peak_amplitudes_ave_ied = max(abs(average_ieds(:, search_start:search_end)), [], 2);
            writematrix(peak_amplitudes_ave_ied, strrep(ied_timing_file, '.txt', '_peak_amplitudes_50ms_ave_IED.txt'), 'Delimiter', 'tab');
       
        end
    end
end

fprintf('Finished processing all subjects.\n');
