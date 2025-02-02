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
        [~, ~, raw_electrode_data] = xlsread(electrode_file_path);

        % Extract headers from raw (first row)
        header_row = raw_electrode_data(1, :);
        
        % Find column indices (allowing variations like x_partial, y_partial, z_partial)
        x_column_index = find(strcmpi(header_row, 'x') | strcmpi(header_row, 'x_partial'));
        y_column_index = find(strcmpi(header_row, 'y') | strcmpi(header_row, 'y_partial'));
        z_column_index = find(strcmpi(header_row, 'z') | strcmpi(header_row, 'z_partial'));

        % Validate indices
        if isempty(x_column_index) || isempty(y_column_index) || isempty(z_column_index)
            fprintf('Skipping %s %s: One or more coordinate columns (X/Y/Z or X_partial/Y_partial/Z_partial) not found.\n', subject_id, run_directories(run_index).name);
            continue;
        end

        % Extract MNI coordinates while preserving N/A values
        mni_coordinates = raw_electrode_data(2:end, [x_column_index, y_column_index, z_column_index]);

        % Convert cell array to numeric while preserving NaN values
        total_contacts = size(mni_coordinates, 1);
        numeric_mni_coordinates = NaN(total_contacts, 3); % Initialize as NaN

        for contact_index = 1:total_contacts
            for coordinate_index = 1:3
                if ischar(mni_coordinates{contact_index, coordinate_index}) % If it's a string
                    if strcmpi(mni_coordinates{contact_index, coordinate_index}, 'N/A') || isempty(mni_coordinates{contact_index, coordinate_index}) % Keep N/A as NaN
                        numeric_mni_coordinates(contact_index, coordinate_index) = NaN;
                    else
                        numeric_mni_coordinates(contact_index, coordinate_index) = str2double(mni_coordinates{contact_index, coordinate_index}); % Convert to number
                    end
                elseif isnumeric(mni_coordinates{contact_index, coordinate_index}) % If already numeric, assign directly
                    numeric_mni_coordinates(contact_index, coordinate_index) = mni_coordinates{contact_index, coordinate_index};
                end
            end
        end

        % Load the transformation matrix as ASCII (since it's not a binary .mat file)
        transformation_matrix = load(transformation_matrix_path, '-ascii');  % Load as numeric matrix

        % Convert to homogeneous coordinates (add a column of ones)
        valid_contacts = ~isnan(numeric_mni_coordinates(:, 1)) & ~isnan(numeric_mni_coordinates(:, 2)) & ~isnan(numeric_mni_coordinates(:, 3));
        homogeneous_mni_coordinates = [numeric_mni_coordinates(valid_contacts, :), ones(sum(valid_contacts), 1)];

        % Apply transformation matrix only for valid numeric rows
        transformed_homogeneous_coordinates = (transformation_matrix * homogeneous_mni_coordinates')';  % Transformed coordinates
        native_space_coordinates = transformed_homogeneous_coordinates(:, 1:3); % Remove the homogeneous coordinate

        % Replace transformed values in the raw data, keeping N/A as is
        valid_contact_counter = 1;
        for contact_index = 1:total_contacts
            if valid_contacts(contact_index)
                raw_electrode_data{contact_index + 1, x_column_index} = native_space_coordinates(valid_contact_counter, 1);
                raw_electrode_data{contact_index + 1, y_column_index} = native_space_coordinates(valid_contact_counter, 2);
                raw_electrode_data{contact_index + 1, z_column_index} = native_space_coordinates(valid_contact_counter, 3);
                valid_contact_counter = valid_contact_counter + 1;
            end
        end


        % Ensure all column names are non-empty and valid
        for column_index = 1:length(header_row)
            if isempty(header_row{column_index}) || ~ischar(header_row{column_index})
                header_row{column_index} = sprintf('Unknown_Column_%d', column_index);
            else
                header_row{column_index} = matlab.lang.makeValidName(header_row{column_index});
            end
        end

        % Convert the updated raw cell array into a table
        output_table = cell2table(raw_electrode_data(2:end, :), 'VariableNames', header_row);

        % Save the table to the new output directory
        output_file_path = fullfile(output_folder, [subject_id '_channel_info_native.xlsx']);
        writetable(output_table, output_file_path, 'FileType', 'spreadsheet');

        % Display completion message
        fprintf('Processed %s %s: Electrode coordinates transformed to native space and saved to %s.\n', subject_id, run_directories(run_index).name, output_file_path);
    end
end

disp('All subjects and runs processed.');
