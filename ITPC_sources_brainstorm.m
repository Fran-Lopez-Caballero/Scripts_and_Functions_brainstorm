%% Compute ITPC in source scouts 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script computes ITPC in every vertex of a given scout in brainstorm
% using iterations in paralel with parfor. Complex spectra for every vertex
% and single trial is extracted from brainstorm and ITPC is calculated from
% it using Brian Coffman's function. Then, ITPCs from all vertices within a
% scout are averaged and resulting values are stored in external directory

% Instructions: 

% 1) Introduce subject codes in "Define subjects" section
% 2) Change conditions, scouts, etc. inside the parfor loop (lines 35 to 50)
% 3) Run script and see results in brainstorm (reload folders from subject)

% Written by Fran López-Caballero for CNRL on 05/17/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define subjects

participant = {'xxxx','xxxx'}; % Subject codes

%% Compute ITPC in signal averaged across scout vertices

tic
disp(' ');      
disp('-------------------------');  
disp('COMPUTING ITPC IN SPECIFIC SOURCES');  
disp(datetime)
disp('-------------------------');
disp(' '); 

parfor p = 1:length(participant)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Define variables inside the loop %%%%%%%%%%%%%%%%%%%%
    root_dir_bs = '/data/CNRL/data04/brainstorm_db/PEPP_Intensity_EAGBR'; % Where your protocol is
    save_dir = '/data/CNRL/data04/ITPC_source_results'; % Directory where output data will be extracted
    protocolName = 'PEPP_Intensity_EAGBR';
    condition = {'75dB_nf_notch','80dB_nf_notch','85dB_nf_notch'}; % Name of condition folders
    average_file_string = 'notch_average_'; % Unique identifier in name of file where average is stored under condition folder in BS
    kernel_name = 'results_MN_MEG'; % Same than before for kernel
    history_trial_identifier = 'data_'; % string identifying trial name in History variable (leave like this)
    Freq_window = 20:80; % In Hz % 20:80
    which_atlas = 'HCPMMP1'; % Atlas where used scouts are ('Intensity_AC_New')
    scout_list = {'L_A1_ROI L','R_A1_ROI R'}; % Extracts a single file with vertices from these scouts
    % Scouts must have exactly the same name as in subject anatomy: load tess_cortex_pial_low.mat from 
    % that subject's anatomy (anat folder) in matllab, and check label names undder Atlas(x).Scouts variable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%% Load brainstorm in each thread %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addpath('/data/CNRL/data03/PEPP/Fran/Global_functions');
    addpath('~/matlab/brainstorm3_v20220706');
    % '/data/CNRL/data03/PEPP/Fran/Global_functions'
    % 'X:/PEPP/Fran/Global_functions'
    brainstorm server
    iProtocol = bst_get('Protocol', protocolName);
    if isempty(iProtocol)
        error(['Unknown protocol: ' protocolName]);
    end
    % Select the current procotol
    gui_brainstorm('SetCurrentProtocol', iProtocol);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load anatomy to know scout vertices later %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Outcome = load([root_dir_bs '/anat/' participant{p} '/brainstormsubject.mat'],'Cortex');
    if isempty(Outcome)
        error([participant{p} ' has no anatomy']);
    end
    Cortex_file = load([root_dir_bs '/anat/' Outcome.Cortex],'Atlas');
    Atlas = Cortex_file.Atlas;
    Pos_main_atlas = find(strcmp({Atlas.Name},which_atlas));
    if isempty(Pos_main_atlas)
        error(['There is no ' which_atlas ' atlas in participant ' participant{p}]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Reload subject first
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);

    % Prepare cell array to store sources ITPC
    ITPC_table = {};

    for c = 1:length(condition)
        % Generate good trial files
        folders = dir([root_dir_bs '/data/' participant{p} '/' condition{c} '/']);
        if isempty(folders) % Some subjects don't have all frequencies
            warning([participant{p} ' has no files for ' condition{c}]);
            continue;
        end
        average_file = find(contains({folders.name},average_file_string));
        if isempty(average_file)
            error([participant{p} ' ' condition{c} ' has no average']);
        elseif length(average_file) > 1
            error([participant{p} ' ' condition{c} ' has more than one average']);
        end
        % Load history from average to retrieve trials composing it
        Outcome = load([root_dir_bs '/data/' participant{p} '/' condition{c} '/' folders(average_file).name],'History');
        History = Outcome.History;
        % Retrieve trial list
        pos_list = find(contains(History(:,3),'List of averaged files'));
        if (isempty(pos_list)) || (length(pos_list) > 1)
            error([participant{p} ': cannot identify trials for average in History file']);
        end
        raw_trials_names = History(pos_list+1:end,3);
        if isempty(raw_trials_names)
            error([participant{p} ': cannot identify trials for average in History file']);
        end
        
        % Get current kernel name
        kernel_file = find(startsWith({folders.name},kernel_name) & ~endsWith({folders.name},'copy.mat'));
        if isempty(kernel_file)
            error([participant{p} ' ' condition{c} ' has no source kernel']);
        elseif length(kernel_file) > 1
            error([participant{p} ' ' condition{c} ' has no source kernel']);
        end
        kernel_bs_name = [participant{p} '/' condition{c} '/' folders(kernel_file).name];

        % Get trial names and combine them with temporal kernel name
        trials = {}; pos = 1;
        for i = 1:length(raw_trials_names)
            if ~contains(raw_trials_names{i},'_trial')
                error([participant{p} ': something odd with trial names']);
            end
            trials{pos} = ['link|' kernel_bs_name '|' participant{p} '/' condition{c} '/' raw_trials_names{i}(strfind(raw_trials_names{i},history_trial_identifier):end)];
            pos = pos + 1;
        end
        % Example output: 
        % 'link|2044C/80dB_nf_notch/results_MN_MEG_GRAD_MEG_MAG_KERNEL_230308_1007.mat|2044C/80dB_nf_notch/data_80dB_nf_notch_trial002_bl.mat'

        % In case some conditions are missing, continue to next condition    
        if isempty(trials)
            continue;
        end

        for sl = 1:length(scout_list)
            
            disp(' ');      
            disp('-------------------------');
            disp(['Extracting single-trial complex spectra for ' participant{p} ' condition ' condition{c} ' scout ' scout_list{sl}]);
            disp(datetime)
            disp(' ');  

            % Process: Time-frequency (Morlet wavelets)
            sFiles_single_trials = bst_process('CallProcess', 'process_timefreq', trials, [], ...
                'clusters',      {which_atlas, {scout_list{sl}}}, ...
                'scoutfunc',     5, ...  % All
                'edit',          struct(...
                     'Comment',         'Scouts,Complex,20-80Hz', ...
                     'TimeBands',       [], ...
                     'Freqs',           Freq_window, ...
                     'MorletFc',        1, ...
                     'MorletFwhmTc',    3, ...
                     'ClusterFuncTime', 'none', ...
                     'Measure',         'none', ...
                     'Output',          'all', ...
                     'SaveKernel',      0), ...
                'normalize2020', 0, ...
                'normalize',     'none');  % None: Save non-standardized time-frequency maps

            % Process: Extract values: 
            sFiles_extracted = bst_process('CallProcess', 'process_extract_values', sFiles_single_trials, [], ...
                'timewindow', [], ...
                'freqrange',  [Freq_window(1), Freq_window(end)], ...
                'rows',       '', ...
                'isabs',      0, ...
                'avgtime',    0, ...
                'avgrow',     0, ...
                'avgfreq',    0, ...
                'matchrows',  1, ...
                'dim',        2, ...  % Concatenate time (dimension 2)
                'Comment',    '');

            % Once the data is extracted, delete ITPC from each trial
            sFiles = bst_process('CallProcess', 'process_delete', sFiles_single_trials, [], ...
                'target', 1);  %#ok<*NASGU> % Delete selected files

            % Calculate ITPC from complex spectra (Brian's function)
            data=load([root_dir_bs '/data/' sFiles_extracted.FileName]);
            S = length(data.TFmask(1,:)); % S = number of samples
            T = length(data.TF(1,:,1))/S; % T = Number of trials
            C = length(data.RowNames); % C = Number of channels/vertices
            F = length(data.Freqs); % F = Number of frequencies
            ComplexSpctrm = permute(data.TF,[3 1 2]); % Change dimorder to F C S*T
            ComplexSpctrm = reshape(ComplexSpctrm,[F C S T]); % Redim to F C S T
            ComplexSpctrm = permute(ComplexSpctrm,[4 2 3 1]); % Change dimorder to T C S F 
            ITPC=data;
            ITPC.Time = ITPC.Time(1:S);

            disp(' ');      
            disp('-------------------------');
            disp(['Computing ITPC for ' participant{p} ' condition ' condition{c} ' scout ' scout_list{sl}]);
            disp(datetime)
            disp(' ');  

            % compute inter-trial phase coherence (itpc)
            ITPC.TF       = ComplexSpctrm./abs(ComplexSpctrm);         % divide by amplitude
            ITPC.TF       = squeeze(abs(mean(ITPC.TF))); % this will give the itc

            % Add comment and measure to see in Brainstorm 
            ITPC.Comment = ['ITPC_' condition{c} '_' scout_list{sl}];
            ITPC.Measure = 'ITPC';

            % Incorporate generated variable to brainstorm
            temp = [participant{p} '/' condition{c} '/channel_vectorview306_acc1.mat']; % grab any file to select folder (channel one, for instance)
            [sStudy, iStudy] = bst_get('AnyFile', temp);    
            OutputFile = db_add(iStudy,ITPC);

            % Delete large file in brainstorm
            if ~isempty(sFiles_extracted.FileName)
                delete([root_dir_bs '/data/' sFiles_extracted.FileName])
            end

            % Average vertices within scout
            sFiles = bst_process('CallProcess', 'process_average_rows', OutputFile, [], ...
                'avgtype',   3, ...  % Average by scout (if signal names include vertex index)
                'avgfunc',   1, ...  % Arithmetic average: mean(x)
                'overwrite', 1);
            
            % Change name of file to retrieve it when extracting ITPC values
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           ['ITPC_' condition{c} '_' scout_list{sl}], ...
                'output',        2);  % Add to file path

            % Delete ITPC from workspace (need to assign to empty in parloop)
            ITPC = {};
            % clear('ITPC');

            % Load resulting TF to store in table
            Output = load([root_dir_bs '/data/' sFiles.FileName],'TF');
            TF = Output.TF;
            % Store information in table, averaged across vertices
            ITPC_table{c+1,1} = [participant{p} '_' condition{c}]; %#ok<*SAGROW>
            if c == 1 % Only first iteration, add header
                ITPC_table{1,sl+1} = scout_list{sl};
            end
            ITPC_table{c+1,sl+1}  = squeeze(TF); % All time points, all frequencies
        end
        
    end
    parsave([save_dir '/' participant{p} '_ITPC_source.mat'],ITPC_table,'ITPC_table');
end

disp 'DONE COMPUTING ITPC IN SPECIFIC SOURCES (RELOAD ALL FOLDERS TO SEE FILES IN BRAINSTORM)';
disp(datetime)
toc
