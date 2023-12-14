% 0) Notes

% This script will prepare anatomy without a FEM
% If willing to get EEG or BIMODAL sources, use Prepare_FEM_anatomy instead
% For this to run, subject_array column 11 must have nothing on it

% Developed for CNRL by Fran López-Caballero
%% 1) Open local brainstorm with the protocol you want to modidy

clear; clc;
% Define protocol
protocol = 'Project_baseline'; % Project_baseline or Project_followup

% Load subject array
root_dir = ['C:/Project/User/' protocol];
load(['C:/Project/User/' protocol '/subject_array']);
participant = {subject_array{:,1}};

% If willing to select an specific participant, uncomment this
% participant = {'xxxx'};

%% 2) Import anatomy folder with HCP

tic
disp(' ');      
disp('-------------------------');  
disp('IMPORTING ANATOMY DATA FOR Project');
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    % Only for those without anatomy already imported
    if ~strcmp(subject_array{pos_subj,11},'needs_atlas') && ~strcmp(subject_array{pos_subj,11},'needs_merge') && ~strcmp(subject_array{pos_subj,11},'Anatomy_ready') 
        try
            
            % Reload anatomy folder of that subject too
            disp(' ');      
            disp('-------------------------');
            disp(['loading anatomy folder for participant ' participant{p}]);
            disp(datetime)
            disp(' ');

            prot_subs = bst_get('ProtocolSubjects');
            current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
            db_reload_subjects(current_sub);
            
            % Input files
            sFiles = [];
            folder = dir(['C:/Project/HCPproc/' participant{p}]); % To define date
            infolder = find(~contains({folder.name},'.') & ~contains({folder.name},'dontuse')); % The ones that are not dots
            if length(infolder) > 1 % Because of the stupid "_do_not_use" files...
                error(['More than one HCP folder file for ' participant{p}]);
            end
            RawFiles = ['C:/Project/HCPproc/' participant{p} '/' folder(infolder).name '/T1w/' participant{p}];
            if ~exist(['C:/Project/HCPproc/' participant{p} '/' folder(infolder).name '/T1w/' participant{p}],'file')
                error(['No anatomy folder for ' participant{p}]);
            end
            % Also, show error if no hcp labels
            dir_labels = dir(['C:/Project/HCPproc/' participant{p} '/' folder(infolder).name '/T1w/' participant{p} '/label/']);
            hcp_labs = find(contains({dir_labels.name},'HCPMMP1'));
            hcp_labs2 = find(contains({dir_labels.name},'hcpmmp1'));
            if isempty(hcp_labs) && isempty(hcp_labs2)
                error(['No HCP labels in ' participant{p} ' HCProc folder']);
            end
            % lh.HCPMMP1.annot rh.HCPMMP1.annot

            % Process: CNRL import anatomy folder
            sFiles = bst_process('CallProcess', 'process_CNRL_import_anatomy', sFiles, [], ...
                'subjectname', participant{p}, ...
                'mrifile',     {RawFiles, 'FreeSurfer'}, ...
                'nvertices',   75000, ... % Changed to 75,000 as in conversation in Teams on 11/07/22 
                'nas',         [0, 0, 0], ...
                'lpa',         [0, 0, 0], ...
                'rpa',         [0, 0, 0], ...
                'ac',          [0, 0, 0], ...
                'pc',          [0, 0, 0], ...
                'ih',          [0, 0, 0], ...
                'aseg',        1); 
            
            % If successful, update subject_array for this subject
            subject_array{pos_subj,11} = 'needs_atlas';
            save([root_dir '/subject_array.mat'],'subject_array')
        catch
            % If unsuccessful, update subject_array for this subject
            subject_array{pos_subj,11} = 'Problems';
            save([root_dir '/subject_array.mat'],'subject_array')
        end
    end
end

disp 'DONE WITH IMPORTING ANATOMY DATA FOR Project!!!'
disp(datetime)
toc

%% 3) Define fiducials and ensure MNI coordinates are set (manual)
pause
%% 4) Create personalized atlas for each subject based on HCPMMP1 one

% Directory specification
anat_dir = ['C:/Project/brainstorm_db/' protocol '/anat/']; 

% If willing to select an specific participant, uncomment this
% participant = {'xxxx'};
% Only if wanting to do so in the common anatomy, uncomment this
% participant = {'@default_subject'};

% Example of atlas 1:
title_new_atlas = 'Auditory_and_frontal'; 
% List of scouts in HCPMMP1 atlas that will be used to create new one
new_scout_list = {'L_52_ROI L','L_A1_ROI L','L_A4_ROI L','L_A5_ROI L','L_LBelt_ROI L',...
    'L_MBelt_ROI L','L_OP4_ROI L','L_PBelt_ROI L','L_RI_ROI L','L_STSdp_ROI L','L_STSvp_ROI L',...
    'R_52_ROI R','R_A1_ROI R','R_A4_ROI R','R_A5_ROI R','R_LBelt_ROI R','R_MBelt_ROI R',...
    'R_OP4_ROI R','R_PBelt_ROI R','R_RI_ROI R','R_STSdp_ROI R','R_STSvp_ROI R', 'L_45_ROI L',...
    'L_IFSa_ROI L','L_IFSp_ROI L','R_45_ROI R','R_IFSa_ROI R','R_IFSp_ROI R'...
    'L_OFC_ROI L','L_pOFC_ROI L','R_OFC_ROI R','R_pOFC_ROI R'};

new_scout_names = {'52_L','A1_L','A4_L','A5_L','LBelt_L', ...
    'MBelt_L','OP4_L','PBelt_L','RI_L','STSdp_L','STSvp_L', ...
    '52_R','A1_R','A4_R','A5_R','LBelt_R','MBelt_R', ...
    'OP4_R','PBelt_R','RI_R','STSdp_R','STSvp_R','45_L', ...
    'IFSa_L','IFSp_L','45_R','IFSa_R','IFSp_R' ...
    'OFC_L','pOFC_L','OFC_R','pOFC_R'};

% In same order, list of new colors of scouts in new atlas
new_scout_colors = {[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]...
    [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]...
    [0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]};

% Now Call function 
[subject_array] = local_function_personalize_scouts(root_dir, anat_dir, participant,title_new_atlas,new_scout_list,new_scout_names,new_scout_colors, subject_array);

%% 5) Merge scouts of interest

% Directory specification
anat_dir = ['C:/Project/brainstorm_db/' protocol '/anat/']; 

% If willing to select an specific participant, uncomment this
% participant = {'xxxx'};
% Only if wanting to do so in the common anatomy, uncomment this
% participant = {'@default_subject'};

which_atlas = 'Auditory_and_frontal'; % Atlas from which scouts will be merged

% As many as needed
scouts_to_merge_group = {{'A1_L','LBelt_L','PBelt_L'},{'A1_R','LBelt_R','PBelt_R'},...
    {'A1_L','LBelt_L','MBelt_L','PBelt_L'},{'A1_R','LBelt_R','MBelt_R','PBelt_R'},...
    {'45_L','IFSa_L','IFSp_L'},{'45_R','IFSa_R','IFSp_R'},...
    {'OFC_L', 'pOFC_L'},{'OFC_R', 'pOFC_R'}}; % Exactly as they appear in Atlas
name_merged_scout_group = {'AUDCORTEX3ROI_L','AUDCORTEX3ROI_R','AUDCORTEX_L','AUDCORTEX_R',...
    'IFG_L','IFG_R','OFC_L_merged','OFC_R_merged'};

merged_scout_color = [0,0,0]; % in RGB
remove_original_scouts = 0; % 0 = NO; 1 = YES; % remove or not the original individual scouts merged

% Now Call function 
[subject_array] = local_function_merge_scouts(root_dir, anat_dir, participant,which_atlas,scouts_to_merge_group,name_merged_scout_group,merged_scout_color,remove_original_scouts, subject_array);

%% Local function stored (don't need to touch)

function [subject_array] = local_function_personalize_scouts(root_dir, anat_dir, participant,title_new_atlas,new_scout_list,new_scout_names,new_scout_colors, subject_array)

for p = 1:length(participant)
    % Sources-or-anatomy-specific problems
    if ~strcmp(participant{p},'@default_subject')
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        if ~strcmp(subject_array{pos_subj,11},'needs_atlas')
            continue; % on to next subject
        end
    else
        disp('Changing default anatomy scouts');
    end
    if ~exist([anat_dir participant{p} '/brainstormsubject.mat'],'file');continue;end
    load([anat_dir participant{p} '/brainstormsubject.mat'])    
    load([anat_dir Cortex]) % contains the atlas
    pos_new = find(strcmp({Atlas.Name},title_new_atlas)); %#ok<*EFIND>
    if isempty(pos_new) % We have to create this atlas
        Pos_main_atlas = find(strcmp({Atlas.Name},'HCPMMP1'));
        % Ojo que puede ser hcpmmp1 también
        if isempty(Pos_main_atlas)
            Pos_main_atlas = find(strcmp({Atlas.Name},'hcpmmp1'));
        end
        if isempty(Pos_main_atlas)
            error(['no atlas under HCPMMP1 or hcpmmp1 names found for ' participant{p}])
        end
        % Crate new entry in Atlas with the name chosen
        Pos_new_atlas = size(Atlas,2)+1;
        Atlas(Pos_new_atlas).Name = title_new_atlas;  %#ok<*AGROW>
        % For every scout selected for this atlas, add it to the new entry
        for nal = 1:length(new_scout_list)
            Pos_orig = find(strcmp({Atlas(Pos_main_atlas).Scouts.Label},new_scout_list{nal}));
            if isempty(Pos_orig)
                error([new_scout_list{nal} ' label not found in original atlas']);
            end
            Atlas(Pos_new_atlas).Scouts(nal) = Atlas(Pos_main_atlas).Scouts(Pos_orig); %#ok<*FNDSB>
            Atlas(Pos_new_atlas).Scouts(nal).Label = new_scout_names{nal};
            Atlas(Pos_new_atlas).Scouts(nal).Color = new_scout_colors{nal}/255;
        end
        % Make it the default atlas
        iAtlas = Pos_new_atlas; %#ok<*NASGU>
        % Identify variables in file to be sure we save back the same
        variableInfo = who('-file',[anat_dir Cortex]);
        save([anat_dir Cortex],variableInfo{:})
        % Identify subject as completed
        if ~strcmp(participant{p},'@default_subject')
            subject_array{pos_subj,11} = 'needs_merge';
            save([root_dir '/subject_array.mat'],'subject_array')
        end
        
    end
end

% Reload subject
for p = 1:length(participant)
    if ~strcmp(participant{p},'@default_subject')
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        % Sources-or-anatomy-specific problems
        if ~strcmp(subject_array{pos_subj,11},'needs_atlas')
            continue; % on to next subject
        end
    else
        disp('Changing default anatomy scouts');
    end
    try
        prot_subs = bst_get('ProtocolSubjects');
        if strcmp(participant{p},'@default_subject') % default anatomy
            current_sub = 0;
        else % Any subject
            current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
        end
        db_reload_subjects(current_sub);
    catch
        disp(['Please reload anatomy for participant ' participant{p} ' to see the new atlas'])
    end
end

end

%% Local function stored (don't need to touch)

function [subject_array] = local_function_merge_scouts(root_dir, anat_dir, participant,which_atlas,scouts_to_merge_group,name_merged_scout_group,merged_scout_color,remove_original_scouts, subject_array)
for p = 1:length(participant)
    % Sources-or-anatomy-specific problems
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if ~strcmp(participant{p},'@default_subject')
        if ~strcmp(subject_array{pos_subj,11},'needs_merge')
            continue; % on to next subject
        end
    else
        disp('Changing default anatomy scouts');
    end
    
    % Once for every set of scouts we want to merge
    for scomerge = 1:length(name_merged_scout_group)
    
        % Define this merged_scout and its components
        name_merged_scout = name_merged_scout_group{scomerge};
        scouts_to_merge = scouts_to_merge_group{scomerge};
        
        if ~exist([anat_dir participant{p} '/brainstormsubject.mat'],'file');continue;end
        load([anat_dir participant{p} '/brainstormsubject.mat'])    
        load([anat_dir Cortex]) % contains the atlas
        Pos_main_atlas = find(strcmp({Atlas.Name},which_atlas));
        if isempty(Pos_main_atlas)
            error(['There is no ' which_atlas ' atlas in participant ' participant{p}]);
        end
        % Find out if the scout you are about to create already exists
        newscout_exists = find(strcmp({Atlas(Pos_main_atlas).Scouts.Label},name_merged_scout), 1);
        if ~isempty(newscout_exists)
            % If it exists, continue to next subject
            continue;
        end
        % Create an empty struct with appropiate fields to store the scouts to merge
        f = fieldnames(Atlas(Pos_main_atlas).Scouts)';
        f{2,1} = {};
        matrix_scouts_to_merge = struct(f{:});   
        % Add the scout information to the newly created struct
        for i = 1:length(scouts_to_merge)
            Pos_scout = find(strcmp({Atlas(Pos_main_atlas).Scouts.Label},scouts_to_merge{i}));
            if isempty(Pos_scout)
                error(['There is no ' scouts_to_merge{i} ' scout in ' which_atlas ' atlas for participant ' participant{p}]);
            end
            matrix_scouts_to_merge(i) = Atlas(Pos_main_atlas).Scouts(Pos_scout);
        end

        % === Join scouts ===
        pos_to_add = length(Atlas(Pos_main_atlas).Scouts)+1;
        Atlas(Pos_main_atlas).Scouts(pos_to_add).Seed = matrix_scouts_to_merge(1).Seed; %#ok<*AGROW>
        Atlas(Pos_main_atlas).Scouts(pos_to_add).Region = matrix_scouts_to_merge(1).Region;
        Atlas(Pos_main_atlas).Scouts(pos_to_add).Handles = [];
        Atlas(Pos_main_atlas).Scouts(pos_to_add).Function = matrix_scouts_to_merge(1).Function;
        try
            Atlas(Pos_main_atlas).Scouts(pos_to_add).Vertices = unique([matrix_scouts_to_merge.Vertices]);
        catch
            % Try again flipping the matrix before
            for i = 1:length(matrix_scouts_to_merge)
                matrix_scouts_to_merge(i).Vertices = matrix_scouts_to_merge(i).Vertices';
            end
            Atlas(Pos_main_atlas).Scouts(pos_to_add).Vertices = unique([matrix_scouts_to_merge.Vertices]);
        end
        Atlas(Pos_main_atlas).Scouts(pos_to_add).Label = name_merged_scout;
        Atlas(Pos_main_atlas).Scouts(pos_to_add).Color = merged_scout_color/255;

        % === Remove old scouts ===
        % RemoveScouts(iScouts);
        if remove_original_scouts == 1
            for i = 1:length(scouts_to_merge)
                Pos_scout = find(strcmp({Atlas(Pos_main_atlas).Scouts.Label},scouts_to_merge{i}));
                Atlas(Pos_main_atlas).Scouts(Pos_scout) = []; 
            end
        end

        % Save back surface file with new Atlas/scouts
        variableInfo = who('-file',[anat_dir Cortex]);
        % Identify variables in file to be sure we save back the same
        save([anat_dir Cortex],variableInfo{:})
    
    end
    
    disp(' ');
    disp('-------------------------');  
    disp(['Scouts merged for ' participant{p} ': reload its anatomy to see new merged scouts']);
    disp(datetime)
    disp(' ');
    
    % Identify subject as completed
    if ~strcmp(participant{p},'@default_subject')
        subject_array{pos_subj,11} = 'Anatomy_ready';
        save([root_dir '/subject_array.mat'],'subject_array')
    end
    
end
end
