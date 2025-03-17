% arc_single_subj_FC_analysis  Performs EEG-fMRI analysis for one subject & run on ARC cluster
% without requiring SPM. Uses MATLAB's built-in niftiread/niftiinfo instead.
%
% (a) Extract the average BOLD time series for each electrode (3×3×3 voxel neighborhood),
% (b) Compute seed-based static or dynamic FC for each main channel,
% (c) Examine the relationship between EEG amplitude and fMRI FC,
% (d) Regress out electrode distance from the seed.
%
% Written by Tahereh Rashnavadi, Feb 2025



function arc_extract_mean_fMRI_TS(subjectID, runID, baseEEGDir, baseFMriDir, baseElecDir)
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
    outDir = fullfile(baseEEGDir, 'seed_staticFC_results', subjectID, runID);
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
    
    %     % If your MATLAB version doesn't handle .nii.gz directly, we can unzip it first:
    %     unzippedFile = gunzip_if_needed(fmriFile);
    
    %% 5. Extract BOLD Time Series for each electrode (3x3x3 neighborhood)
    % Extracts fMRI timeseries using FSL's fslmeants within MATLAB
    % Uses a (3×3×3)-1 voxel neighborhood around each electrode location
    % Averages the 26 voxels to mitigate signal loss due to electrode-induced artifacts
    % Ensures robust BOLD signals for functional connectivity (FC) analysis
    
    % --- % DEFINE OUTPUT FILE FOR EXTRACTED TIMESERIES FOR THIS SUBJECT/RUN ---
    outDir_TS = fullfile(baseEEGDir, 'fMRI_timeseries', subjectID, runID);
    if ~exist(outDir_TS, 'dir')
        mkdir(outDir_TS);
    end
    
    electrodeTS = nan(nVolumes, nDepthElectrodes); % Correctly preallocate for time series
    
    % Define coordinates file for fslmeants (zero-based indexing for FSL)
    % Compute 3×3×3 neighborhood for each electrode (zero-based for FSL)
    for i = 1:nDepthElectrodes
        % Extract the center (zero-based) coordinates for this electrode
        x = depthCoords(i,1);
        y = depthCoords(i,2);
        z = depthCoords(i,3);
        % IMPORTANT: This assumes your nativeCoords are already zero-based. If your coordinates are 1-based, then you either need to subtract 1 again or adjust accordingly.
        % -------------------------------------------------------------------------
        timeseriesFile = fullfile(outDir_TS, sprintf('%s_%s_fMRI_ts.txt', subjectID, runID));
        maskFile = fullfile(outDir_TS, sprintf('roi_%d.nii.gz', i));
    
        % Step 1: Create a 3×3×3 Neighborhood Mask
        fslCmd1 = sprintf('fslmaths %s -roi %d 3 %d 3 %d 3 0 -1 %s', fmriFile, x, y, z, maskFile);
        system(fslCmd1);
    
        % Step 2: Extract Mean Time Series from the Mask
        fslCmd2 = sprintf('fslmeants -i %s -m %s -o %s', fmriFile, maskFile, timeseriesFile);
        [status, cmdout] = system(fslCmd2);
    
        if status ~= 0
            warning('fslmeants failed for electrode %s at (%d, %d, %d): %s', depthChannelNames(i), x, y, z, cmdout);
            continue;
        end
    
        % Step 3: Load Extracted Time Series
        tsData = load(timeseriesFile);
        fprintf('Extracted time series for electrode %s at (%d, %d, %d):\n', depthChannelNames(i), x, y, z);
    
        % Ensure time series length matches expected fMRI volume count
        if length(tsData) ~= nVolumes
            warning('Mismatch in time series length for electrode %s. Expected %d, got %d.', ...
                depthChannelNames(i), nVolumes, length(tsData));
            continue;
        end
    
        electrodeTS(:, i) = tsData;
        % Remove mask file after use**
        if exist(maskFile, 'file')
            delete(maskFile);
        end
    end
    % Define the file paths
    matFilePath = fullfile(outDir_TS, sprintf('%s_%s_electrodeTS.mat', subjectID, runID));
    csvFilePath = fullfile(outDir_TS, sprintf('%s_%s_electrodeTS.csv', subjectID, runID));
    
    % Save as .MAT file (MATLAB format)
    save(matFilePath, 'electrodeTS', 'depthChannelNames');
    
    % Convert electrodeTS matrix to a table with headers
    electrodeTS_table = array2table(electrodeTS, 'VariableNames', depthChannelNames);
    
    % Save as CSV with headers
    writetable(electrodeTS_table, csvFilePath);
    
   
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

function outFile = gunzip_if_needed(gzFile)
% If your MATLAB supports reading .nii.gz directly, you may skip this step
% and pass gzFile to niftiread. Otherwise, we unzip it if needed.
    [filePath, fileName, fileExt] = fileparts(gzFile);
    if strcmpi(fileExt, '.gz')
        outFileCandidate = fullfile(filePath, fileName); % .nii
        if ~exist(outFileCandidate,'file')
            fprintf('Unzipping %s...\n', gzFile);
            gunzip(gzFile, filePath);
        end
        outFile = outFileCandidate;
    else
        outFile = gzFile; % already a .nii or no .gz extension
    end
end


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
