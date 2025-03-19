% arc_single_subj_FC_analysis  Performs EEG-fMRI analysis for one subject & run on ARC cluster
% without requiring SPM. Uses MATLAB's built-in niftiread/niftiinfo instead.
%
% (a) Extract the average BOLD time series for each electrode (3×3×3 voxel neighborhood),
% (b) Compute seed-based static or dynamic FC for each main channel,
% (c) Examine the relationship between EEG amplitude and fMRI FC,
% (d) Regress out electrode distance from the seed.
%
% Written by Tahereh Rashnavadi, Feb 2025



function arc_single_subj_analysis_all_electrodes(subjectID, runID, iedTypes, baseEEGDir, baseFMriDir, baseElecDir)
    %ARC_SINGLE_SUBJ_FC_ANALYSIS  Performs EEG-fMRI analysis for one subject & run.
    %
    % Input arguments:
    %   subjectID   (char) - e.g. 'ICE001'
    %   runID       (char) - e.g. 'Run1'
    %   iedTypes    (cell) - e.g. {'IED1','IED2'}
    %   baseEEGDir  (char) - base directory for EEG data
    %   baseFMriDir (char) - base directory for fMRI data
    %   baseElecDir (char) - base directory for electrode coords
    %
    
    % Example usage: ON MY HARD DRIVE
    % addpath('/Users/trashnavadi/Documents/2024/postdoc/Levan/analysis/arc_cluster/scripts');
    % arc_single_subj_FC_analysis('ICE001', 'Run1', ...
    %     {'IED1','IED2'}, ...                    % iedTypes
    %     '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara', % baseEEGDir: where the Peak EEG of Averaged IEDs of each channel are saved
    %     '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/ICE_denoised_filtered_funcs', ... % baseFMriDir
    %     '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/native_space_electrodes_coords' % baseElecDir: in Native space as the fMRI images
    % );
    
    % Example usage, ARC CLUSTER
    % addpath('/work/levan_lab/Tara')
    % arc_single_subj_FC_analysis('ICE001', 'Run1', ...
    %     {'IED1','IED2'}, ...                    % iedTypes
    %     '/work/levan_lab/Tara', ...             % baseEEGDir
    %     '/work/goodyear_lab/Tara/ICE_denoised_filtered_funcs', ... % baseFMriDir
    %     '/work/levan_lab/Tara/native_space_electrodes_coords' ...  % baseElecDir
    % );
    
    %% Setup: (No longer adding SPM path)
    % If you need FSL or something else, do so here, e.g. setenv('FSLDIR', '...')
    % addpath(genpath('/usr/local/fsl/etc/matlab'))  % Example for FSL's MATLAB code if needed
    %% 1. make an output directory
    % --- CREATE AN OUTPUT DIRECTORY FOR THIS SUBJECT/RUN ---
    outDir = fullfile(baseEEGDir, 'all_electrodes_as_seeds_staticFC_results', subjectID, runID);
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    %% 2. Load electrode coordinates (native space) from the excel sheet
    electrodeFile = fullfile(baseElecDir, [subjectID '_native_space_coords.xlsx']);
    if ~exist(electrodeFile,'file')
        warning('Electrode file not found for %s. Skipping.', subjectID);
        return;
    end
    
    T = readtable(electrodeFile);
    requiredCols = {'ElectrodeName_labConvention_', 'ContactNumber', 'x', 'y', 'z', 'ElectrodeType'};
    
    if ~all(ismember(requiredCols, T.Properties.VariableNames))
        error('Expected columns not found in %s', electrodeFile);
    end
    
    % Use the actual column header for electrode names
    channelNames = string(T.ElectrodeName_labConvention_) + string(T.ContactNumber); % Create full channel name = e.g. 'dLA' + '1' => 'dLA1'
    nativeCoords = round([T.x, T.y, T.z]); % Round to nearest voxel
    electrodeTypes = T.("ElectrodeType");   % Column specifying if electrode is depth or grid/strip
    
    % Remove Electrodes with Missing Coordinates (NaN)
    validIdx = all(~isnan(nativeCoords), 2);
    if any(~validIdx)
        missingChannels = channelNames(~validIdx);
        warning('The following electrodes are missing coordinates and will be excluded: %s', strjoin(missingChannels, ', '));
    end
    
    % Keep only valid electrodes
    channelNames = channelNames(validIdx);
    nativeCoords = nativeCoords(validIdx, :);
    electrodeTypes = electrodeTypes(validIdx);
    
    nElectrodes    = size(nativeCoords,1); % # of depth electrodes
    
    %% 3. Identify Depth Electrodes
    depthElectrodeIdx = find(strcmpi(electrodeTypes, 'depth'));  % Get indices of depth electrodes
    if isempty(depthElectrodeIdx)
        warning('No depth electrodes found for subject %s.', subjectID);
        return;
    end
    
    depthChannelNames = channelNames(depthElectrodeIdx);  % Names of depth electrodes
    depthCoords       = nativeCoords(depthElectrodeIdx, :);  % Coordinates of depth electrodes
    nDepthElectrodes = size(depthCoords, 1);
    
    %% 4. Load fMRI data (.nii or .nii.gz) with built-in MATLAB I/O
    fmriFile = fullfile(baseFMriDir, subjectID, runID,[subjectID '_denoised_filtered_func_' runID '.nii.gz']);
    
    if ~exist(fmriFile,'file')
        warning('fMRI file not found: %s. Skipping.', fmriFile);
        return;
    end
    
    % Use `fslinfo` to get details about the fMRI file
    fslInfoCmd = sprintf('fslinfo %s', fmriFile);
    [status, infoOutput] = system(fslInfoCmd);
    
    if status ~= 0
        warning('Failed to retrieve fMRI info for %s. Check FSL installation.', fmriFile);
        return;
    end
    
    % Parse `fslinfo` output
    infoLines = strsplit(infoOutput, '\n'); % Split output into lines
    nVolumes = [];
    TR = [];
    
    for ll = 1:length(infoLines)
        line = strtrim(infoLines{ll});
        tokens = regexp(line, '\s+', 'split'); % Split by tabs or multiple spaces
    
        % Look for 'dim4' and extract its numeric value
        if length(tokens) >= 2 && strcmp(tokens{1}, 'dim4')
            nVolumes = str2double(tokens{2});
        elseif length(tokens) >= 2 && strcmp(tokens{1}, 'pixdim4')
            TR = str2double(tokens{2});
        end
    end
    
    % Ensure values were extracted correctly
    if isempty(nVolumes) || isnan(nVolumes)
        error('Could not determine number of volumes (dim4) from fMRI file: %s', fmriFile);
    end
    if isempty(TR) || isnan(TR)
        error('Could not determine TR (pixdim4) from fMRI file: %s', fmriFile);
    end
    
    fprintf('fMRI file: %s\n - Number of volumes: %d\n - TR: %.2f sec\n', fmriFile, nVolumes, TR);
       
    %% 5. Call the previously extracted BOLD Time Series for each electrode (3x3x3 neighborhood)
    % RUN arc_extract_mean_fMRI_TS.m to Extracts fMRI timeseries using FSL's fslmeants within MATLAB
    % Uses a (3×3×3)-1 voxel neighborhood around each electrode location
    % Averages the 26 voxels to mitigate signal loss due to electrode-induced artifacts
    % Ensures robust BOLD signals for functional connectivity (FC) analysis
    
    % --- % DEFINE OUTPUT FILE FOR EXTRACTED TIMESERIES FOR THIS SUBJECT/RUN ---
    mean_fMRI_TS_path = fullfile(baseEEGDir, 'fMRI_timeseries', subjectID, runID, sprintf('%s_%s_electrodeTS.mat', subjectID, runID));

    % Check if the pre-extracted file exists
    if ~exist(mean_fMRI_TS_path, 'file')
        error('Pre-extracted fMRI time series file not found: %s', mean_fMRI_TS_path);
    end

    % Load the pre-extracted electrode time series
    load(mean_fMRI_TS_path, 'electrodeTS', 'depthChannelNames');

    % Verify loaded data
    if isempty(electrodeTS) || isempty(depthChannelNames)
        error('Loaded fMRI time series data is empty or incomplete: %s', mean_fMRI_TS_path);
    end

    fprintf('Loaded pre-extracted fMRI time series for %s, Run: %s\n', subjectID, runID);
    fprintf('Number of electrodes: %d, fMRI Volumes: %d\n', length(depthChannelNames), size(electrodeTS,1));

    % **Recreate the electrodeTS_table from the loaded data**
    electrodeTS_table = array2table(electrodeTS, 'VariableNames', depthChannelNames);

    % (Optional) Save the recreated table as CSV if needed
    csvFilePath = fullfile(baseEEGDir, 'fMRI_timeseries', subjectID, runID, sprintf('%s_%s_electrodeTS.csv', subjectID, runID));
    writetable(electrodeTS_table, csvFilePath);
   
     %% 6. Compute Static FC (Pearson correlation)
    % Find electrodes with all-zero time series
    zeroElectrodes = all(electrodeTS == 0, 1);
    
    % Check if any electrodes have all-zero values
    if any(zeroElectrodes)
        removedChannels = depthChannelNames(zeroElectrodes);
    
        % Explicitly print the removed electrodes
        fprintf('The following electrode contacts are out of the functional image (all-zero time series) and will be removed:\n');
        disp(removedChannels);  % Display as a cell array
    
        % Ensure the warning message is displayed
        warning('Removing %d electrodes: %s', length(removedChannels), strjoin(removedChannels, ', '));
    end
    
    % Remove zero-value electrodes from analysis
    electrodeTS(:, zeroElectrodes) = [];
    depthChannelNames(zeroElectrodes) = []; % Also update depth electrode names
    nDepthElectrodes = size(electrodeTS, 2); % Update electrode count
    
    % --- staticFC matrix  ---
    staticFC = corr(electrodeTS, 'Rows', 'complete');  % Compute Pearson correlation
   
    % --- SAVE the staticFC matrix to .mat ---
    staticFCFile = fullfile(outDir, sprintf('%s_%s_staticFC.mat', subjectID, runID)); % Saves the static functional connectivity matrix to a .mat file in outDir.
    save(staticFCFile, 'staticFC');
    
    % Plot and Save FC Matrix as Image
    figure;
    colormap jet;
    imagesc(staticFC);
    colorbar;
    title(sprintf('Static FC (z-transformed) Matrix - %s %s', subjectID, runID));
    xlabel('Electrodes');
    ylabel('Electrodes');
    set(gca, 'XTick', 1:nDepthElectrodes, 'XTickLabel', depthChannelNames, 'XTickLabelRotation', 45);
    set(gca, 'YTick', 1:nDepthElectrodes, 'YTickLabel', depthChannelNames);
    grid on;
    
    % Save FC Image
    fcImgFile = fullfile(outDir, sprintf('%s_%s_staticFC.png', subjectID, runID));
    saveas(gcf, fcImgFile);
    
    fprintf('Static FC matrix saved as: %s\n', fcImgFile);
    
    %% 7. Load Peak EEG amplitudes of IED average per channel (for an example IED type)
    if isempty(iedTypes)
        warning('No iedTypes specified. Skipping amplitude loading.');
        return;
    end

    % Define new filenames for saving results
    fcAnalysisMatFile = fullfile(outDir, sprintf('%s_%s_FC_analysis_results.mat', subjectID, runID));
    fcAnalysisExcelFile = fullfile(outDir, sprintf('%s_%s_FC_analysis_results.xlsx', subjectID, runID));

    % Initialize result storage
    results = struct();
    data_table = {};
    
    % Initialize a structure to store EEG amplitudes for all IED types
    EEG_amplitudes_by_IED = struct();
    
    for t = 1:length(iedTypes)
        iedType = iedTypes{t}; % Get current IED type

        % Construct the EEG amplitude file path for this IED type
        IED_eegAmplitudeFile = fullfile(baseEEGDir, subjectID, 'peak_amp_IED_per_channel', ...
            [subjectID '_' runID '_' iedType '_adjusted_peak_amplitudes_50ms_ave_IED.txt']);

        if ~exist(IED_eegAmplitudeFile, 'file')
            warning('EEG amplitude file not found for %s, run %s, IED type %s. Skipping.', subjectID, runID, iedType);
            continue;
        end

        % From EEG LOG: Read EEG amplitude table for this specific IED type
        IED_ave_info = readtable(IED_eegAmplitudeFile, 'Delimiter', 'tab');

        % Ensure required columns are present
        if ~ismember('Channel', IED_ave_info.Properties.VariableNames)
            warning('Column ''Channel'' not found in EEG amplitude file for subject %s.', subjectID);
            continue;
        end
        if ~ismember('PeakAmplitude', IED_ave_info.Properties.VariableNames)
            error('Column ''PeakAmplitude'' not found in %s.', IED_eegAmplitudeFile);
        end

        % Extract EEG channel names and amplitudes
        channelNamesEEG = string(IED_ave_info.Channel);

        [commonChannels, fMRI_Idx, EEG_Idx] = intersect(depthChannelNames, channelNamesEEG);
        fprintf('Total common electrodes between EEG and fMRI for %s: %d\n', iedType, length(commonChannels));
        
        % Debugging: Check if final indices match
        fprintf('Final matched electrodes: %d (should match staticFC size: %d and EEG size: %d)\n', ...
            length(fMRI_Idx), size(EEG_Idx, 1), length(commonChannels));
        
        % Ensure the indices match the corresponding EEG table
        yIED = IED_ave_info.PeakAmplitude(EEG_Idx); % Extract EEG amplitudes of averaged IEDs per channel for matched electrodes

        % Store results in the struct using IED type as the key
        EEG_amplitudes_by_IED.(iedType) = yIED;   
        
        % Define the key for this subject and IED type
        thisKey = sprintf('%s_%s', subjectID, iedType);

        % Get all the EEG channels for this IED type
        if isempty(commonChannels)
            warning('No common channels found for subject %s, IED type %s. Skipping.', subjectID, iedType);
            continue;
        end
    
        % Assign all commonChannels instead of just the main ones
        mainChNames = commonChannels;

        % Find indices of these main channels in commonChannels
        [tf, mainIdx] = ismember(mainChNames, commonChannels);

        fprintf('EEG amplitude of averaged-IED for all channels stored for %s with %d channels.\n', iedType, length(mainChNames)); 

        mainIdx = mainIdx(tf); % Keep only valid ones
          
        fprintf('Processing IED Type: %s, Subject: %s\n', iedType, subjectID);
        disp(mainChNames);
    
        %% 9. Perform Seed-Based FC Analysis for Main Channels
        % EEG channles have already been updated by filtering out the
        % channels that were depth electrodes but did not have coords in
        % the excel sheet--> matchedEEGChannels
        % Find common electrodes between matchedEEGChannels and fMRI/depth electrode channels, 
        % the number of common electrodes could be less than each of fMRI or EEG
        % channesl due to some channels not have fMRI timeseires (voxels out of image bound) which might be recorded by EEG
        % and this makes the number of common channels less than EEG
        % channels too, also there are channels that EEG has been recorded
        % there but their coordinates are not provided in the excel sheet.
        
        for m = 1:length(mainIdx)
            seedIdx = mainIdx(m);
            seedName = commonChannels{seedIdx};
    
            if ~strcmp(seedName, depthChannelNames(fMRI_Idx(seedIdx)))
                warning('Seed electrode %s (Idx=%d) is not valid for fMRI analysis since the channel contacts were out of the brain functional boundary and removed.', seedName, seedIdx);
                continue;
            end
    
            fprintf('\n--- Seed: %s (Idx=%d), IED Type: %s ---\n', seedName, seedIdx, iedType);
    
            % (i) Extract Seed Time Series
            % Extract Seed Time Series by Matching the Column Header Name
            seedTS = electrodeTS_table{:, seedName};  % Extract seed time series as a numeric column vector
    
            % (ii) Extract Time Series of Other Depth Electrodes (ONLY those in commonChannels)
            commonChannelsTS = electrodeTS_table{:, commonChannels};  % Extract numeric matrix for FC calculation
    
            % (iii) Compute Functional Connectivity (FC)
            fcValues = corr(seedTS, commonChannelsTS);  % Compute Pearson correlation
            fcValues(strcmp(commonChannels, seedName)) = 1;  % Self-correlation (optional)
            fcValues_Z = atanh(fcValues); % Fisher’s Z-transformation
    
            % (iv) Compute Distance from Seed to common Depth Electrodes
            seedIdx = find(strcmp(commonChannels, seedName));  % Find index for distance calculation
            % Step 1: Find indices of survived common channels in depthChannelNames to get their coordinates
            [~, commonDepthIdx] = ismember(commonChannels, depthChannelNames);  % Indices of surviving channels
            commonDepthCoords = depthCoords(commonDepthIdx, :);  % Extract coordinates of the survived electrodes
    
            % Step 2: Find the coordinates of the seed electrode
            seedCoord = nativeCoords(strcmp(depthChannelNames, seedName), :);  % Get seed's coordinates
    
            % Step 3: Compute Euclidean distance from seed to all surviving electrodes
            distVals = sqrt(sum((commonDepthCoords - seedCoord).^2, 2));  % Distance computation

             % ** Exclude the seed electrode itself (where FC = 1) to avoid
             % bias, also it avoids atanh(1) which is inf to be in
             % fcvalues_Z
            validIdx = fcValues ~= 1;  % Mask to remove seed electrode itself as it acts as an outlier
            % Remove the seed electrode from all variables 
            xFC_filtered = fcValues_Z(validIdx);
            xFC_filtered = xFC_filtered(:);
            yIED_filtered = yIED(validIdx);
            dist2seed_filtered = distVals(validIdx);

            % Exclude electrodes that are too close (<= threshold) and 
            % (v) Compare FC with EEG Amplitude, Controlling for Distance
            % A minimum of four data points is required to compute partial correlation reliably.
%             threshold_dist = 6;
            threshold_dist = 7;
            validIdx_dist = dist2seed_filtered >= threshold_dist;

            % Check if there are any valid electrodes left
            if ~any(validIdx_dist)
                warning('All electrodes are closer to the seed than the threshold (%d mm). Stopping analysis.', threshold_dist);
                continue; % Stops execution of the current function/script
            end

            xFC_valid_distance = xFC_filtered(validIdx_dist); % FC values to depth electrodes
            yIED_valid_distance = yIED_filtered(validIdx_dist); % EEG values to depth electrodes
            dist2seed_valid = dist2seed_filtered(validIdx_dist);
     
            
            % **Check if there are at least 4 data points for correlation and regression**
            num_valid_electrodes = numel(xFC_valid_distance);
            if num_valid_electrodes < 4
                warning('Not enough valid electrodes (%d) after filtering. Skipping this seed.', num_valid_electrodes);
                continue; % Skip correlation and regression
            end
            
            % Pearson correlation: FC vs. EEG amplitude (Without Distance Correction)
            r_simple   = corr(xFC_valid_distance, yIED_valid_distance);

            % Partial Correlation (Correcting for Distance)
            r_partial  = partialcorr(xFC_valid_distance, yIED_valid_distance, dist2seed_valid);

            % **Apply Fisher z-transformation to correlation values**
            z_simple   = atanh(r_simple);
            z_partial  = atanh(r_partial);
    
            % Store correlation results
            results(t).seeds(m).name = seedName;
            results(t).seeds(m).IED_type = iedType;
            results(t).seeds(m).r_simple = r_simple;
            results(t).seeds(m).r_partial = r_partial;
            results(t).seeds(m).z_simple = z_simple;  % Fisher z-transformed correlation
            results(t).seeds(m).z_partial = z_partial;
   
            % (vi) Multiple Regression: IED ~ FC + logDist
            % Dependent Variable: EEG amplitude (IED). Independent Variables: FC and log-transformed distance (logDist).
            logDist_filtered = log(dist2seed_valid + eps); % Add eps to avoid log(0)
            %         invDist_filtered = 1 ./ dist2seed_valid; % Avoids division by zero since distances are > 0

            % **Ensure valid regression data before creating the table**
            if isempty(xFC_valid_distance) || isempty(logDist_filtered) || isempty(yIED_valid_distance)
                warning('Regression variables are empty after filtering. Skipping this seed.');
                continue;
            end

            % Create a regression table with the interaction term
            tbl = table(xFC_valid_distance, logDist_filtered, yIED_valid_distance, 'VariableNames', {'FC','LogDist', 'IED'});
    
            % Fit the regression model with the interaction term
            mdl = fitlm(tbl, 'IED ~ FC + LogDist');
    
            % **Extract Beta Coefficients from the Model**
            beta0 = mdl.Coefficients.Estimate(1); % Intercept
            beta1 = mdl.Coefficients.Estimate(2); % FC effect
            beta2 = mdl.Coefficients.Estimate(3); % Log(Distance) effect
            % Extract p-values for each coefficient (Intercept, FC, LogDist)
            p_value_intercept = mdl.Coefficients.pValue(1); % p-value for Intercept
            p_value_fc = mdl.Coefficients.pValue(2);        % p-value for FC
            p_value_logDist = mdl.Coefficients.pValue(3);   % p-value for Log(Distance)

    
            % Store results in structured format
            results(t).seeds(m).beta0 = beta0;
            results(t).seeds(m).beta1 = beta1;
            results(t).seeds(m).beta2 = beta2;
    
            % Compute ANOVA table for extracting F-statistic and p-value
            anova_results = anova(mdl);
    
            % Store R², RMSE, and model performance metrics
            results(t).seeds(m).R_squared = mdl.Rsquared.Ordinary;
            results(t).seeds(m).Adjusted_R_squared = mdl.Rsquared.Adjusted;
            results(t).seeds(m).RMSE = mdl.RMSE;
    
            % Extract F-statistic and p-value from ANOVA table
            results(t).seeds(m).F_stat = anova_results.F(2); % Extracts F-statistic for the model
            results(t).seeds(m).p_value = anova_results.pValue(2); % Extracts p-value for the model
            results(t).seeds(m).p_value_intercept = p_value_intercept;
            results(t).seeds(m).p_value_fc = p_value_fc;
            results(t).seeds(m).p_value_logDist = p_value_logDist;

            F_stat = anova_results.F(2); % Extract F-statistic
            p_value = anova_results.pValue(2); % Extract p-valu
    
            % Append to table format for Excel, Extract Akaike Information Criterion (AIC) for model comparison
            newRow = {subjectID, runID, iedType, seedName, r_simple, z_simple, r_partial, z_partial, ...
                beta0, beta1, beta2, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted, mdl.RMSE, ...
                F_stat, p_value, p_value_fc, p_value_logDist, mdl.ModelCriterion.AIC};

            data_table = [data_table; newRow];
    
            % **Display Regression Model Results**
            fprintf('Seed %s: simple corr(FC, IED) = %.3f, partial corr = %.3f\n', seedName, z_simple, z_partial);
            fprintf('\nRegression for seed %s:\n', seedName);
            disp(mdl);
    
            %% ** Visualizations**
            % make matlab not to display the figures
            set(0, 'DefaultFigureVisible', 'off');

            cleanKey = strrep(thisKey, '_', ' '); % Replace underscore with space or ''
    
            % (a) Scatter Plot: FC vs. EEG Amplitude wrt distance**
            % Normalize distance for better visualization of size and color
            minSize = 70; % Minimum marker size
            maxSize = 200; % Maximum marker size
            % Close electrodes = Big markers, Far electrodes = Small markers
            markerSizes = maxSize - (maxSize - minSize) * (dist2seed_valid - min(dist2seed_valid)) / (max(dist2seed_valid) - min(dist2seed_valid));
    
            % Define colormap (e.g., "jet" colormap for smooth transition)
            colormap jet;
            c = dist2seed_valid; % Color based on distance
    
            % Create scatter plot with color and size encoding distance
            figure('Name', sprintf('%s, FC vs EEG Amplitude - %s', cleanKey, seedName));
            scatter(xFC_valid_distance, yIED_valid_distance, markerSizes, c, 'filled', 'MarkerEdgeColor', 'k'); % Size & color map to distance
            colorbar; % Show color legend for distance
            xlabel('Functional Connectivity (FC)');
            ylabel('IED Amplitude');
            title(sprintf('%s, %s, Seed: %s: FC vs. EEG Amplitude (Color & Size ~ Distance)', cleanKey, runID, seedName));
            grid on;
            box on;
    
            % Save scatter plot
            scatterFile1 = fullfile(outDir, sprintf('%s_%s_%s_FC_vs_EEG_DistanceCoded.png', subjectID, iedType, seedName));
            saveas(gcf, scatterFile1);
    
            % **(b) Scatter Plot: FC vs. Distance**
            figure('Name', sprintf('%s, Seed:  %s, FC vs Distance', cleanKey, seedName));
            scatter(dist2seed_valid, xFC_valid_distance, 80, 'b', 'filled', 'MarkerEdgeColor', [0 0 0.5]);
            xlabel('Distance to Seed Electrode');
            ylabel('Functional Connectivity (FC)');
            title(sprintf('%s, %s, Seed %s: FC vs. Distance', cleanKey, runID, seedName));
            grid on;
            box on;

            scatterFile2 = fullfile(outDir, sprintf('%s_%s_%s_FC_vs_Distance.png', subjectID, iedType, seedName));
            saveas(gcf, scatterFile2);
    
            % **(c) Scatter Plot: EEG Amplitude of IEDs vs. Distance**
            figure('Name', sprintf('%s, seed: %s, Distance_vs_IED', cleanKey, seedName));
            scatter(dist2seed_valid, yIED_valid_distance, 100, [1 0.5 0], 'filled', 'MarkerEdgeColor', [0.8 0.3 0]);
            xlabel('Distance to Seed Electrode');
            ylabel('IED Amplitude');
            title(sprintf('%s, %s, Seed: %s, EEG Amplitude of IEDs vs. Distance', cleanKey, runID, seedName));
            grid on;
            box on;

            scatterFile3 = fullfile(outDir, sprintf('%s_%s_%s_Distance_vs_IED.png', subjectID, iedType, seedName));
            saveas(gcf, scatterFile3);
    
            % Best for: Simple, clear comparisons of FC across electrodes.
%             figure('Name', sprintf('%s, Seed: %s, Seed-Based FC Line Plot', cleanKey, seedName));
%             stem(fcValues, 'filled');
%             xticks(1:length(commonChannels));
%             xticklabels(commonChannels);
%             xtickangle(45);
%             xlabel('Electrodes');
%             ylabel('Functional Connectivity');
%             title(sprintf('%s, %s, Seed: %s, Seed-Based FC', cleanKey, runID,seedName));
%             grid on;
        end
    end

    % Save the results as a .mat file
    save(fcAnalysisMatFile, 'results');

    % Save the results as an Excel file

    % Define column names for results table
    column_names = {'Subject', 'Run', 'IED Type', 'Seed', 'Simple Corr (FC, IED)', ...
    'Z_transformed Simple Corr (FC, IED)', 'Partial Corr (FC, IED)', ...
    'Z_transformed Partial Corr (FC, IED)', 'Beta0 (Intercept)', 'Beta1 (FC)', ...
    'Beta2 (LogDist)', 'R-Squared', 'Adjusted R-Squared', 'RMSE', ...
    'F-Statistic', 'p-value (Overall Model)', 'p-value (FC)', 'p-value (LogDist)', 'AIC'};

   
    % Initialize data_table properly before appending rows
    if ~exist('data_table', 'var') || isempty(data_table)
        data_table = cell(0, numel(column_names)); % Ensure it has the right structure
    end

    if ~isempty(data_table)
        results_table = cell2table(data_table, 'VariableNames', column_names);
    else
        warning('No valid data points were processed. Table will not be created.');
        results_table = table(); % Return an empty table to avoid errors
    end

    writetable(results_table, fcAnalysisExcelFile);

    fprintf('Results saved for subject %s, run %s in %s\n', subjectID, runID, outDir);


    %% ----------------------------------------------------------------
    % (v) Dynamic FC (Optional)
    % For demonstration: 20-sec window
    %             windowSizeSec = 90;
    %             windowSizeVol = round(windowSizeSec / TR);
    %             stepSize      = 1; % shift by 1 volume
    %
    %             dynFC = compute_dynamic_fc(electrodeTS, seedIdx, windowSizeVol, stepSize);
    % ----------------------------------------------------------------


      
%% ----------------------------------------------------------------
% 10.  HELPER FUNCTIONS
%% ----------------------------------------------------------------

% function outFile = gunzip_if_needed(gzFile)
% % If your MATLAB supports reading .nii.gz directly, you may skip this step
% % and pass gzFile to niftiread. Otherwise, we unzip it if needed.
%     [filePath, fileName, fileExt] = fileparts(gzFile);
%     if strcmpi(fileExt, '.gz')
%         outFileCandidate = fullfile(filePath, fileName); % .nii
%         if ~exist(outFileCandidate,'file')
%             fprintf('Unzipping %s...\n', gzFile);
%             gunzip(gzFile, filePath);
%         end
%         outFile = outFileCandidate;
%     else
%         outFile = gzFile; % already a .nii or no .gz extension
%     end
% end


% function dynamicFC = compute_dynamic_fc(electrodeTS, seedIndex, windowSize, stepSize)
% %COMPUTE_DYNAMIC_FC Computes seed-based dynamic connectivity with a sliding window.
% %
% %   dynamicFC = compute_dynamic_fc(electrodeTS, seedIndex, windowSize, stepSize)
% % 
% % Inputs:
% %   electrodeTS - [#Volumes x #Electrodes]
% %   seedIndex   - which electrode index to treat as seed
% %   windowSize  - how many volumes per window
% %   stepSize    - shift in volumes for each successive window
% %
% % Output:
% %   dynamicFC   - [#Windows x #Electrodes] matrix of correlation values
% 
%     [nVolumes, nElectrodes] = size(electrodeTS);
%     nWindows = floor((nVolumes - windowSize)/stepSize) + 1;
%     dynamicFC = zeros(nWindows, nElectrodes);
%     
%     seedDataFull = electrodeTS(:, seedIndex);
%     
%     for w = 1:nWindows
%         startVol = (w - 1)*stepSize + 1;
%         endVol   = startVol + windowSize - 1;
%         
%         seedWindow = seedDataFull(startVol:endVol);
%         for elecIndex = 1:nElectrodes
%             otherWindow = electrodeTS(startVol:endVol, elecIndex);
%             dynamicFC(w, elecIndex) = corr(seedWindow, otherWindow, 'type', 'Pearson');
%         end
%     end
% end

end
