%% Transforming the electrodes coordinates from MNI space into native space (ICE study)
% written by Tahereh Rashnavadi, Feb 2025

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
for subject_index = 1:length(subject_list)
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
        transformation_matrix_path = fullfile(run_folder_path, 'standard2highres.mat');
        
        if ~isfile(transformation_matrix_path)
            fprintf('Skipping %s %s: No transformation matrix found.\n', subject_id, run_directories(run_index).name);
            continue;
        end
        
        % Load the Excel file (raw to keep text and numbers)
        [~, ~, raw_data] = xlsread(electrode_file_path);

        % Extract headers from raw (first row)
        header_row = raw_data(1, :);
        
        % Find column indices for X, Y, Z coordinates
        x_idx = find(strcmpi(header_row, 'x'));
        y_idx = find(strcmpi(header_row, 'y'));
        z_idx = find(strcmpi(header_row, 'z'));

        % Ensure coordinate columns are found
        if isempty(x_idx) || isempty(y_idx) || isempty(z_idx)
            fprintf('Skipping %s %s: X, Y, or Z coordinate columns missing.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        % Extract all data rows (excluding headers)
        data_rows = raw_data(2:end, :);  

        % Extract MNI coordinates
        mni_coords = data_rows(:, [x_idx, y_idx, z_idx]);

        % Convert to numeric format, preserving NaN values
        num_contacts = size(mni_coords, 1);
        numeric_mni_coords = NaN(num_contacts, 3); % Initialize as NaN

        for i = 1:num_contacts
            for j = 1:3
                if ischar(mni_coords{i, j}) % If it's a string
                    if strcmpi(mni_coords{i, j}, 'N/A') || isempty(mni_coords{i, j})
                        numeric_mni_coords(i, j) = NaN; % Keep N/A as NaN
                    else
                        numeric_mni_coords(i, j) = str2double(mni_coords{i, j}); % Convert to number
                    end
                elseif isnumeric(mni_coords{i, j}) % If already numeric
                    numeric_mni_coords(i, j) = mni_coords{i, j};
                end
            end
        end

        % Load transformation matrix
        T = load(transformation_matrix_path, '-ascii');  

        % Convert to homogeneous coordinates (add a column of ones)
        valid_rows = ~isnan(numeric_mni_coords(:, 1)) & ~isnan(numeric_mni_coords(:, 2)) & ~isnan(numeric_mni_coords(:, 3));
        homogeneous_coords = [numeric_mni_coords(valid_rows, :), ones(sum(valid_rows), 1)];

        % Apply transformation matrix only to valid coordinates
        transformed_coords_h = (T * homogeneous_coords')';  
        native_coords = transformed_coords_h(:, 1:3); % Remove homogeneous coordinate

        % Replace transformed values in the raw data while preserving all other columns
        counter = 1;
        for i = 1:num_contacts
            if valid_rows(i)
                data_rows{i, x_idx} = native_coords(counter, 1);
                data_rows{i, y_idx} = native_coords(counter, 2);
                data_rows{i, z_idx} = native_coords(counter, 3);
                counter = counter + 1;
            end
        end

        % Ensure the column "Electrode Name (lab convention)" is preserved
        electrode_name_idx = find(strcmpi(header_row, 'Electrode Name (lab convention)'));
        if isempty(electrode_name_idx)
            fprintf('Warning: Electrode Name column missing for %s %s.\n', subject_id, run_directories(run_index).name);
        end

        % Ensure valid column names for the table
        for col_idx = 1:length(header_row)
            if isempty(header_row{col_idx}) || ~ischar(header_row{col_idx})
                header_row{col_idx} = sprintf('Unknown_Column_%d', col_idx);
            else
                header_row{col_idx} = matlab.lang.makeValidName(header_row{col_idx});
            end
        end

        % Convert updated raw data into a table, preserving all columns
        output_table = cell2table(data_rows, 'VariableNames', header_row);

        % Save the transformed data
        output_file_path = fullfile(output_folder, [subject_id '_channel_info_native.xlsx']);
        writetable(output_table, output_file_path, 'FileType', 'spreadsheet');

        % Print completion message
        fprintf('Processed %s %s: Electrode coordinates transformed and saved to %s.\n', subject_id, run_directories(run_index).name, output_file_path);
    end
end

disp('All subjects and runs processed.');
