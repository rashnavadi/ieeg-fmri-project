%% Transforming the electrodes coordinates from MNI space into native space (ICE study)
% Reads depth electrode coordinates in MNI space from Excel.
% Transforms them into the subject’s native anatomical space (T1).
% Further transforms them into the subject’s functional fMRI space.
% Saves the final coordinates as an Excel file (.xlsx).

% written by Tahereh Rashnavadi, Feb 2025

% ON MY HARD DRIVE 
% ICE_reg_folder = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs';
% electrodes_coords_folder  = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/original_ICE/coordinates';
% output_folder = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/native_space_electrodes_coords'; % baseElecDir: in Native space as the fMRI images

% Define base directories
ICE_reg_folder = '/work/goodyear_lab/Tara/ICE_denoised_filtered_funcs/'; % Path to subject reg folders
electrodes_coords_folder = '/work/levan_lab/eegfmri_epilepsy/coordinates/'; % Path to original electrode coordinate files
output_folder = '/work/levan_lab/Tara/native_space_electrodes_coords/'; % New directory for transformed coordinates

% Ensure output folder exists
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% List of subjects (ICE001 to ICE070)
subject_list = arrayfun(@(subject_number) sprintf('ICE%03d', subject_number), 1:70, 'UniformOutput', false);

% Loop through each subject
for subject_index = 13:length(subject_list) % from subject 13 i am analysing
    subject_id = subject_list{subject_index};
    
    % Find corresponding electrode coordinates file
    electrode_file_path = fullfile(electrodes_coords_folder, [subject_id '_channel_info.xlsx']);
    
    if ~isfile(electrode_file_path)
        fprintf('Skipping %s: No electrode coordinates file found.\n', subject_id);
        continue;
    end
    
    % Get list of runs for this subject
    subject_folder_path = fullfile(ICE_reg_folder, subject_id);
    run_directories = dir(fullfile(subject_folder_path, 'Run*')); % Find all runs (Run1, Run2a, Run3b, etc.)

    % Loop through each run
    for run_index = 1:length(run_directories)
        run_folder_path = fullfile(subject_folder_path, run_directories(run_index).name, 'reg');
        % Load transformation matrix from MNI to native anatomical space
        mni2anat_matrix_path = fullfile(run_folder_path, 'standard2highres.mat');
        
        if ~isfile(mni2anat_matrix_path)
            fprintf('Skipping %s %s: No MNI to native transformation matrix found.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        % Load the transformation matrix
        mni2anat_mat = load(mni2anat_matrix_path, '-ascii'); % 4x4 affine transformation matrix

        % Load the Excel file (raw to keep text and numbers)
        [~, ~, raw_data] = xlsread(electrode_file_path);

        % Extract headers from raw (first row)
        header_row = raw_data(1, :);

        % All data rows (exclude header)
        all_rows = raw_data(2:end, :);  % keep everything
        
        % Find column indices for X, Y, Z coordinates
        x_idx = find(strcmpi(header_row, 'x'));
        y_idx = find(strcmpi(header_row, 'y'));
        z_idx = find(strcmpi(header_row, 'z'));
        type_idx = find(strcmpi(header_row, 'Electrode type')); % Column for depth/non-depth classification


        % Ensure coordinate columns are found
        if isempty(x_idx) || isempty(y_idx) || isempty(z_idx)
            fprintf('Skipping %s %s: X, Y, or Z coordinate columns missing.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        % Create a mask for depth electrodes
        is_depth = strcmpi(all_rows(:, type_idx), 'depth');

        % Extract the subset of rows that are depth electrodes
        depth_rows = all_rows(is_depth, :);

        % If there are no depth electrodes, skip or just save the original data
        if isempty(depth_rows)
            fprintf('No depth electrodes found for %s %s.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        % ====== Convert Depth Rows' X/Y/Z to Numeric ======
        % Extract X, Y, Z data (as a cell array)
        xyz_cells = depth_rows(:, [x_idx, y_idx, z_idx]);

        % Pre-allocate numeric matrix
        num_contacts = size(xyz_cells, 1);
        mni_coords = NaN(num_contacts, 3);

        % Convert each cell to numeric, preserving NaN
        for i = 1:num_contacts
            for j = 1:3
                val = xyz_cells{i, j};
                if ischar(val)
                    % If it's a string, try converting to numeric
                    tmp = str2double(val);
                    if isnan(tmp)
                        % If string is 'N/A' or something not convertible, remain NaN
                        mni_coords(i, j) = NaN;
                    else
                        mni_coords(i, j) = tmp;
                    end
                elseif isnumeric(val)
                    mni_coords(i, j) = val;
                else
                    % If it's an unexpected type (like logical), keep as NaN
                    mni_coords(i, j) = NaN;
                end
            end
        end

        % ====== MNI → Native Space Transformation ======
        % Convert MNI coordinates to homogeneous coordinates (Nx4 matrix)
        homogeneous_mni_coords = [mni_coords'; ones(1, size(mni_coords, 1))]; % Ensure 4xN

        % Apply MNI to native space transformation
        native_coords = (mni2anat_mat * homogeneous_mni_coords)'; % Now multiplication is valid
        native_coords = native_coords(:, 1:3); % Extract only X, Y, Z

        % ====== Check for Native→fMRI Transformation ======
        % Convert from native space to fMRI space
        anat2func_matrix_path = fullfile(run_folder_path, 'highres2example_func.mat');

        if isfile(anat2func_matrix_path)
            anat2func_mat = load(anat2func_matrix_path,'-ascii');
            homogeneous_native_coords = [native_coords, ones(num_contacts, 1)]';
            fmri_coords = anat2func_mat * homogeneous_native_coords;
            fmri_coords = fmri_coords(1:3, :)'; % Extract only X, Y, Z
        else
            fprintf('Warning: No anatomical to fMRI transformation matrix found. Using native space coordinates.\n');
            fmri_coords = native_coords; % If no fMRI transform, default to native space
        end

        % ====== Put Transformed Coordinates BACK into all_rows ======
        % Convert fmri_coords (Nx3 double) into cell format
        fmri_coords_cell = num2cell(fmri_coords);

        % Overwrite the X/Y/Z columns only for depth electrodes
        all_rows(is_depth, x_idx) = fmri_coords_cell(:, 1);  % updated X
        all_rows(is_depth, y_idx) = fmri_coords_cell(:, 2);  % updated Y
        all_rows(is_depth, z_idx) = fmri_coords_cell(:, 3);  % updated Z

        % ====== Combine headers and all_rows for final output ======
        final_data = [header_row; all_rows];
        
        % Save as a spreadsheet
        transformed_file_path = fullfile(output_folder, sprintf('%s_transformed_coords.xlsx', subject_id));

        % Attempt to convert final_data to a table
        try
            % Using 'VariableNamingRule','preserve' to keep original column names
            output_table = cell2table(final_data(2:end,:), ...
                'VariableNames', string(final_data(1,:)), ...
                'VariableNamingRule', 'preserve');
            
            % Write table to Excel
            writetable(output_table, transformed_file_path, 'FileType','spreadsheet');
        catch
            % Fallback: older MATLAB might not have 'VariableNamingRule'
            writecell(final_data, transformed_file_path);
        end

        fprintf('Saved transformed coordinates (all columns) for %s in %s\n',...
            subject_id, transformed_file_path);
    end
end

disp('All subjects and runs processed.');













