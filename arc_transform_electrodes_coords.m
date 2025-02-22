%% Convert Electrode Coordinates from MNI to Native Space using FSL
% This script:
% - Reads depth electrode coordinates in MNI space from an Excel file.
% - Uses FSL's `img2imgcoord` to transform them into the subject’s native anatomical space.
% - Saves the final coordinates in an Excel file (.xlsx).
%
% Written by: Tahereh Rashnavadi, Feb 2025

%% Define base directories
% ICE_reg_folder = '/work/goodyear_lab/Tara/ICE_denoised_filtered_funcs/'; % Subject registration folders
% electrodes_coords_folder = '/work/levan_lab/eegfmri_epilepsy/coordinates/'; % Original electrode coordinate files
% output_folder = '/work/levan_lab/Tara/native_space_electrodes_coords/'; % Output directory for transformed coordinates

% ON MY HARD DRIVE 
ICE_reg_folder = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs';
electrodes_coords_folder  = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/original_ICE/coordinates';
output_folder = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/native_space_electrodes_coords'; % baseElecDir: in Native space as the fMRI images


% Ensure output folder exists
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% List of subjects (ICE001 to ICE070)
subject_list = arrayfun(@(subject_number) sprintf('ICE%03d', subject_number), 1:70, 'UniformOutput', false);

%% Loop through each subject
for subject_index = 13:length(subject_list)  % Start analysis from ICE013 onwards
    subject_id = subject_list{subject_index};
    
    % Find corresponding electrode coordinates file
    electrode_file_path = fullfile(electrodes_coords_folder, [subject_id '_channel_info.xlsx']);
    
    if ~isfile(electrode_file_path)
        fprintf('Skipping %s: No electrode coordinates file found.\n', subject_id);
        continue;
    end
    
    % Get list of runs for this subject
    subject_folder_path = fullfile(ICE_reg_folder, subject_id);
    run_directories = dir(fullfile(subject_folder_path, 'Run*')); % Find all run folders (Run1, Run2a, Run3b, etc.)

    %% Loop through each run
    for run_index = 1:length(run_directories)
        run_folder_path = fullfile(subject_folder_path, run_directories(run_index).name, 'reg');
        
        % Locate required FSL transformation matrix
        standard2func_matrix_path = fullfile(run_folder_path, 'standard2example_func.mat');
        example_func_path = fullfile(run_folder_path, 'example_func.nii.gz');

        if ~isfile(standard2func_matrix_path)
            fprintf('Skipping %s %s: No standard2example_func.mat found.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        %% Load Electrode Data
        T = readtable(electrode_file_path);

        % Check for required columns
        requiredCols = {'x', 'y', 'z', 'ElectrodeType'};
        if ~all(ismember(requiredCols, T.Properties.VariableNames))
            fprintf('Skipping %s %s: Missing required columns in electrode file.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        % Convert x, y, z columns to numeric (replace 'N/A' with NaN if they are not numeric already)
        if ~isnumeric(T.x)  % Check if x is not numeric
            T.x = str2double(T.x);  % Convert x to numeric
        end
        if ~isnumeric(T.y)  % Check if y is not numeric
            T.y = str2double(T.y);  % Convert y to numeric
        end
        if ~isnumeric(T.z)  % Check if z is not numeric
            T.z = str2double(T.z);  % Convert z to numeric
        end

        % Handle rows that could not be converted to numeric (replace NaN with original data)
        T.x(isnan(T.x)) = NaN;  % Ensure invalid coordinates are NaN
        T.y(isnan(T.y)) = NaN;
        T.z(isnan(T.z)) = NaN;

        % Extract depth electrode coordinates
        is_depth = strcmpi(T.ElectrodeType, 'depth');
        if sum(is_depth) == 0
            fprintf('No depth electrodes found for %s %s.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        % Convert depth electrode coordinates to numeric matrix
        mni_coords = T{is_depth, {'x', 'y', 'z'}};

        % Save MNI coordinates to a text file (.txt) for FSL
        mni_coords_file = fullfile(output_folder, sprintf('%s_MNI_coords.txt', subject_id));
        writematrix(mni_coords, mni_coords_file, 'Delimiter', 'space');

        % Define output file for transformed native-space coordinates
        native_coords_file = fullfile(output_folder, sprintf('%s_native_coords.txt', subject_id));

        %% Run FSL Transformation: MNI → Native Space
        fsl_cmd = sprintf(['/Applications/fsl/share/fsl/bin/img2imgcoord -src /Applications/fsl/data/standard/MNI152_T1_2mm.nii.gz ' ...
            '-dest %s -xfm %s -mm < %s > %s'], ...
            example_func_path, standard2func_matrix_path, mni_coords_file, native_coords_file);

        % Print the FSL command and run it
        fprintf('Running FSL command: %s\n', fsl_cmd);
        [status, cmdout] = system(fsl_cmd);

        %% Load Transformed Native Coordinates
        native_coords = readmatrix(native_coords_file);

        % Ensure native_coords is a matrix with 3 columns (x, y, z)
        if size(native_coords, 2) == 3
            % Check if 'ElectrodeType' contains values that represent depth electrodes (e.g., 'depth' flag)
            depth_rows = strcmp(T.ElectrodeType, 'depth');

            % Create a table with native_coords
            coords_table = array2table(native_coords, 'VariableNames', {'x', 'y', 'z'});

            % Make sure the size of depth_rows matches the number of rows in native_coords
            if sum(depth_rows) ~= size(native_coords, 1)
                fprintf('❌ Error: Number of depth electrodes does not match the number of native coordinates\n');
                return;
            end

            % Assign the new coordinates to the filtered rows of T
            T{depth_rows, {'x', 'y', 'z'}} = table2array(coords_table);

            % Now, preserve original coordinates for non-depth electrodes
            non_depth_rows = ~depth_rows;  % Get the rows where ElectrodeType is not 'depth'

            % Loop through the non-depth electrodes and keep their original coordinates
            for i = 1:height(T)
                if non_depth_rows(i)
                    % Preserve original values for x, y, and z for non-depth electrodes
                    T{i, {'x', 'y', 'z'}} = T{i, {'x', 'y', 'z'}};
                end
            end

            % Save updated coordinates as an Excel file
            transformed_file_path = fullfile(output_folder, sprintf('%s_transformed_coords.xlsx', subject_id));
            writetable(T, transformed_file_path, 'FileType', 'spreadsheet');

            fprintf('✅ Saved transformed coordinates for %s in %s\n', subject_id, transformed_file_path);
            
            %% Remove the intermediate files after saving transformed coordinates
            delete(mni_coords_file);  % Delete MNI coordinates file
            delete(native_coords_file);  % Delete native coordinates file
        else
            fprintf('❌ Error: Transformed coordinates do not have three columns (x, y, z)\n');
        end
    end
end


disp('🎯 All subjects and runs processed successfully.');


