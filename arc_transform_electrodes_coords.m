%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Batch std2imgcoord from MNI → Native Space for multiple subjects
% 
% This script:
%  1) Loops over a list of subjects.
%  2) For each subject, reads an Excel file containing coordinates in MNI space.
%  3) Extracts x,y,z columns and writes them to a temporary text file.
%  4) Calls `std2imgcoord` to transform MNI coords → native space (example_func).
%  5) Reads back the transformed coordinates, places them into the table,
%     and writes out a new Excel file with updated coordinates.
%
% Written by: Tahereh Rashnavadi, Feb 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Define directories and subject list
ICE_reg_folder = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs';
input_coords_folder = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/original_ICE/coordinates';
output_coords_folder = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/native_space_electrodes_coords';

% Make sure the output folder exists
if ~exist(output_coords_folder, 'dir')
    mkdir(output_coords_folder);
end

% Define subjects - for example, ICE001 through ICE070
subject_list = arrayfun(@(n) sprintf('ICE%03d', n), 1:70, 'UniformOutput', false);

%% 2) Loop over subjects
for iSubj = 13:length(subject_list)
    
    subject_id = subject_list{iSubj};
    fprintf('\n=== Processing %s ===\n', subject_id);
    
    % Path to subject's Excel file (with MNI coords)
    excel_file = fullfile(input_coords_folder, [subject_id '_channel_info.xlsx']);
    
    % Skip if the file doesn't exist
    if ~isfile(excel_file)
        fprintf('  --> No Excel file found for %s. Skipping.\n', subject_id);
        continue;
    end

    % Read the full Excel table
    T = readtable(excel_file);

    % Check for required columns
    requiredCols = {'x', 'y', 'z', 'ElectrodeType'};
    if ~all(ismember(requiredCols, T.Properties.VariableNames))
        fprintf('Skipping %s: Missing required columns in electrode file.\n', subject_id);
        continue;
    end

    %% Convert x, y, z columns to numeric if they aren't already
    % (e.g., if they're stored as strings in Excel)
    if ~isnumeric(T.x)
        T.x = str2double(T.x);
    end
    if ~isnumeric(T.y)
        T.y = str2double(T.y);
    end
    if ~isnumeric(T.z)
        T.z = str2double(T.z);
    end

    % Handle any non-numeric entries that couldn't be converted
    T.x(isnan(T.x)) = NaN;
    T.y(isnan(T.y)) = NaN;
    T.z(isnan(T.z)) = NaN;

    %% Identify depth electrodes
    is_depth = strcmpi(T.ElectrodeType, 'depth');
    if sum(is_depth) == 0
        fprintf('No depth electrodes found for %s. Skipping.\n', subject_id);
        continue;
    end

    % Extract only the MNI coords for depth electrodes
    mni_coords = T{is_depth, {'x', 'y', 'z'}};

    % Create a temporary file for writing MNI coords (for FSL)
    mni_coords_file = fullfile(output_coords_folder, sprintf('%s_MNI_coords.txt', subject_id));
    writematrix(mni_coords, mni_coords_file, 'Delimiter', 'space');

    % Define output file for transformed native-space coordinates
    native_coords_file = fullfile(output_coords_folder, sprintf('%s_native_coords.txt', subject_id));

    
    %----------------------------------------------------------------------
    % 3) Identify the subject's registration files
    %    (Modify the run folder name if needed: e.g., 'Run1a' or 'Run1')
    %----------------------------------------------------------------------
    run_folder = fullfile(ICE_reg_folder, subject_id, 'Run1a', 'reg');
    
    example_func_nii = fullfile(run_folder, 'example_func.nii.gz');
    standard_nii     = fullfile(run_folder, 'standard.nii.gz');
    xfm_mat          = fullfile(run_folder, 'example_func2standard.mat'); % we use the reverse .mat file for std2imgcoord
    
    % Verify needed files exist
    if ~exist(example_func_nii,'file') || ~exist(standard_nii,'file') || ~exist(xfm_mat,'file')
        fprintf('  --> Missing image or matrix for %s. Skipping.\n', subject_id);
        continue;
    end
    
    %----------------------------------------------------------------------
    % 4) Build std2imgcoord command:
    %    We supply standard(=MNI), image(=example_func), and the matrix 
    %    that goes from example_func→standard or standard→example_func
    %----------------------------------------------------------------------
    
    % Output text file (transformed coords in native space)
    tmp_native_txt = fullfile(output_coords_folder, [subject_id '_native_tmp.txt']);
    
    fsl_cmd = sprintf([ ...
        'std2imgcoord ' ...
        '-img "%s" ' ...  % subject image
        '-std "%s" ' ...  % MNI or standard image
        '-xfm "%s" ' ...  % transform matrix
        '%s > %s'], ...   % input coords -> output coords
        example_func_nii, ...
        standard_nii, ...
        xfm_mat, ...
        mni_coords_file, ...
        tmp_native_txt);
    
    fprintf('  --> Running FSL: %s\n', fsl_cmd);
    
    [status, cmdout] = system(fsl_cmd);
    if status ~= 0
        fprintf('  --> std2imgcoord failed for %s:\n%s\n', subject_id, cmdout);
        continue;
    end
    
    %----------------------------------------------------------------------
    % 5) Read back the transformed coords and put them in the table
    %----------------------------------------------------------------------
    if ~isfile(tmp_native_txt)
        fprintf('  --> No transformed file output for %s. Skipping.\n', subject_id);
        continue;
    end
    
    native_coords = readmatrix(tmp_native_txt);
    
    if size(native_coords, 2) ~= 3
        fprintf('  --> Transformed coords do not have 3 columns for %s. Skipping.\n', subject_id);
        continue;
    end
    
    % Place the new (native) coordinates back into T for the depth electrodes only
    T{is_depth, {'x', 'y', 'z'}} = native_coords;
    
    %----------------------------------------------------------------------
    % 6) Save a new Excel file with the updated (native) coords
    %----------------------------------------------------------------------
    out_excel = fullfile(output_coords_folder, [subject_id '_native_space_coords.xlsx']);
    writetable(T, out_excel, 'FileType', 'spreadsheet');
    
    fprintf('  --> Saved transformed coords for %s: %s\n', subject_id, out_excel);
    
    % Clean up temporary files
    delete(mni_coords_file);
    delete(tmp_native_txt);
end

fprintf('\nAll done!\n');
