%% Notes

% Corregistration step, described in data analysis, 
% introduces a transformation in space in the position 
% of the sensors, in order to fit their positions to the 
% head surface we build from the MRI. 

% Sometimes you may want that the transformation you apply 
% to the channel file in one run  or dataset can be applied 
% to the channel file of another run. Normally, brainstorm 
% will copy over the transformation across all channel files 
% within your study that share the same original EEG and MEG 
% channel positions. But, what if the MEG sensor positions 
% are slightly different from one run to another but you 
% still want to use the same transformation? 
% You can use this script
% Keep in mind that what you are copying are "the movements" you make during corregistration from one file to another, not the channel positions per se. 


% Developed for CNRL by Fran López-Caballero

%% Define variables for LLR

root_dir = '/private/USER/private02/HuMon';
root_dir_bs = '/home/cnrl/brainstorm_db/HuMon_LLRseg'; 

% get protocol name, in case we don't run it with server mode
pos_last = find(root_dir_bs == '/', 1, 'last');
ProtocolName = root_dir_bs(pos_last+1:end); clear pos_last

participant = {
'xx'
};
participant_general_list = participant; % so that it has a different name
% Load info of sessions and blocks per subject
load([root_dir '/Brainstorm_pipelines/session_block_array_corregistration.mat'])
load([root_dir '/Brainstorm_pipelines/additional_bad_chans.mat'])
session = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J'}; % default one, will be redefined
condition = {'11' '12' '13' '31' '32' '33' '51' '52' '53' '71' '72' '73' '91' '92' '93'};
shortest_ISI = {'11' '12' '13'};
loudest_condition = {'13' '33' '53' '73' '93'}; % only for stimulus artefact correction
Exp_cond = {'Quietest', 'Medium_dB', 'Loudest', 'Fastest', 'Fast',...
                'Medium_ISI', 'Slow', 'Slowest'};
block = {'b1' 'b2' 'b3'}; % The default one, will actually be adapted with the block_and_runs variable
runs = {'1' '2' '3' '4' '5' '6' '7' '8' '9'}; % same than before
modality_data = {'EEG','MEG','both_mod'};
wave = {'LLR'}; % inherited from previous structure, there is loop with only one value
epoch_wave = {[-0.15, 0.4]}; % LLR only
epoch_wave_bas = {[-0.15, 0.6]}; % LLR, only to create templates for correction of shortest ISI
epoch_baseline = {[-0.05, 0]}; % LLR baseline correction
epoch_covar = {[-3.5, 0.015]}; % from which we calculate noise cov matrices
time_noise_cov = [-3.4, -0.01]; % time window from which to obtain noise covariance values
reject_EEG = [0, 100]; % peak to peak threshold (change name maybe)
reject_EEG_absolute = [0, 50]; % absolute threshold
reject_MEG_grad = [0, 5000];
reject_MEG_mag = [0, 5000];
MLR_highpass = 0; % should be 10, but  already high passsed at preprocessing
MLR_lowpass = 200;
LLR_highpass = 0; % None (because it's already filtered at  0.5 from before)
LLR_lowpass = 20;
tranband_MLR = 3.7; % Transition bandwith, useful to avoid filter artifacts and reduce filter order
tranband_LLR = 0; % Leave at 0 for now
resample_LLR = 500; % LLR resample from 1500 to 500Hz
wave_name = {'low'}; % since both MLR and LLR will only be low_pass in Brainstorm
reref_option = 1; % 0 = NO 1 = YES Yes or no rereference
ref_EEG = 'M1, M2'; % in case we use it, but we won't for now: Alternative: 'AVERAGE'
detrend_option = 2; % 0 = NO detrend; 1 = Brainstorm all epoch; 2 = BVA from first to last 50ms
detrend_option_templates = 0; % 0 = NO detrend; 1 = Brainstorm all epoch; 2 = BVA from first to last 50ms
subtraction_approach =  1; % 0 = No subtraction produced;   1 = YES.
covar_epochs = 1; % 0 = NO covar epochs; 1 = YES, make them
sensor_analysis = 4; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only; 4 = ALL
source_analysis = 2; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only. Cannot be more than one at a time
source_inverse_measure = 'amplitude'; % 'amplitude' or 'dspm2018'
reset_channel_files = 1; % 0 = leave current channel files, 1 = reset them to backup channel files
delete_previous_file = 1; % 1 = delete previous head models and source kernels if reran these steps
ERROR_AVERAGE = {}; % variable to store errors when no average is possible

% Variables for subtraction approach
hann_wind_subtract = 0; % 0 = NO; 1 = YES Apply hanning window to vector to subtract
hann_wind_template = 1; % 0 = NO; 1 = YES Apply half hanning to end of template
hann_start = 0.4; % For now apply from 400 to the end (600 ms). NOT CONSIDERING BASELINE!
all_val = 2; % 1 = all files 2 = all files that are going to be analized
% Note: if you ever change the length of the baseline portion for MLR or
% LLR, change it in the script to correct with subtraction
template = 2; % 1 = block average; 2 = session average; 3 = total average 
% PROBABLY BEST OPTION FOR TEMPLATE IS 2 (Tradeoff between cleanliness of
% templates and affecting MEG data)
delimiter = ',';
startRow = 12;
formatSpec = '%*q%q%q%*s%*s%*s%[^\n\r]';

initialVars = who; % variables up until here, which won't be deleted afterwards
initialVars = who; % twice so that InitialVars itself is not deleted

%% Apply tranformation matrices to other blocks and project EEG to surface in ALL blocks

% IMPORTANT: THIS 'NEW' VERSION IS TO BE USED WHEN MANUAL TRANSFORMATION
% WAS USED ONLY IN FIRST BLOCK FIRST CONDITION (THAT IS, WE SAID 'NO' WHEN
% BRAINSTORM ASKS TO COPY THE TRANSFORMATION TO OTHER CONDITIONS OF THE
% SAME BLOCK, BECAUSE IT FAILS SOMETIMES OR COPIES IT TO OTHER BLOCKS
% SOMETIMES AND SOMETIMES DON'T). SO, FOR CONSISTENCY, ONLY FIRST BLOCK
% FIRST CONDITION IS MODIFIED AND COPIED FORM THERE ANYWHERE ELSE

% IMPORTANT: MODIFIED SO THAT IT WILL ONLY COPY TRANSFORMATIONS IF "MANUAL
% TRANSFORMATION" IS FOUND IN LABELS OF MEG. ALSO, IT WILL COPY ALL MANUAL
% (NOT OTHERS) TRANSFORMATIONS FOUND

% IMPORTANT: IT WILL ALSO PROJECT TO SURFACE ALL BLOCKS (SO IF FIRST BLOCK
% HAS A TRANSFORMATION, IT WILL COPY THAT TO OTHER BLOCKS AND PROJECT EEG TO
% SURFACE IN ALL BLOCKS; IF FIRST BLOCK HAS NO TRANSFORMATION, THEN IT
% WILL ONLY PROJECT EEG SURFACE TO ALL BLOCKS AND FINISH)

tic
disp(' ');      
disp('-------------------------');  
disp('APLYING TRANSFORMATION MATRICES TO OTHER BLOCKS WITHIN SESSION (LLR)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

no_manual_transf = {}; no_transf_pos = 1;

for p = 1:length(participant)
    
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
      
    % Define sessions within this participant
    pos_par = find(strcmp(session_block_array(:,1), participant{p}(1:4)));
    session = session_block_array{pos_par,2}(1,:);

    for s = 1:length(session)
        
        % Define blocks within this session
        pos_ses = find(strcmp(session_block_array{pos_par,2}(1,:), session{s}));
        block = session_block_array{pos_par,2}{2,pos_ses}(1,:);
        
        % Find first block of the session and check if it exist
        if ~exist([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{1} '_normal'], 'dir') % If folder for condition 11 does not exists
            % Crucially, this is the ONLY folder where you are gonna apply
            % a transformation (exceptions would have been noted)
            error([participant{p} '/LLR_11' session{s} '_' block{1} '_normal does not exist'])
        end

        % Open file from first block to obtain transformation matrix
        % If the folder exist this would exist, if not, an error will pop
        % out so no worries
        load([root_dir_bs '/data/' participant{p} '/LLR_11' session{s} '_' block{1} '_normal/channel_vectorview306_acc1.mat'])

        % Find which columns have manual correction (if any)
        [x,y]= find(contains(TransfMegLabels,'manual correction'));
    
        if isempty(x)
            % If no manual transformation is found, do not copy anything to other folders
            % But keep record of it
            no_manual_transf{no_transf_pos,1} = [participant{p} '/LLR_11' session{s} '_' block{1}]; %#ok<*SAGROW>
            no_transf_pos = no_transf_pos +1;

            % And anyway I need to project EEG electrodes to surface!!

            % Search for all other blocks within the session that are not the first
            % (12, 13 and so on from first block are already changed with the gui)
            % Only needed in 'normal' folders: others don't need to be adjusted
            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            for w = 1:length(wave) % always gonna be LLR now
                for c = 1:length(condition) % it's going to have to be applied to all 15 condition folders
                    for b = 1:length(block) % ALSO FIRST BLOCK NEEDS EEG PROJECTION TO SURFACE
                        results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        infolder = find(results);
                        if isempty(infolder) % if no coincidence exists
                            continue 
                        end
                        if size(infolder,2)> 1 % in case of more than one coincidence (should not happen)
                           error(['more than one folder for ' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        end

                        % Get name of any brainstorm file in the folder (from which
                        % channel file will be retrieved by function)
                        dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                        sweeps_norm = contains({dir_sweeps.name},'_trial');
                        sweeps_list = find(sweeps_norm);
                        position = sweeps_list(1);  % first trial found
                        file_name = [participant{p} '/' folders(infolder).name '/' dir_sweeps(position).name];
                        sFiles = file_name;

                        % Project electrodes to surface ONLY
                        sFiles = bst_process('CallProcess', 'process_channel_project', sFiles, []);
                    end
                end
            end

        else % there is at least one manual transformation

            % First, since original block has a transformation, copy to all
            % other blocks within the session

            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            for w = 1:length(wave) % always gonna be LLR now
                for c = 1:length(condition) % it's going to have to be applied to all 15 condition folders
                    for b = 1:length(block) 
                        % All blocks now (because when using the gui we
                        % won't apply to other files, given that freaking
                        % brainstorm sometimes copies them to other blocks
                        % aside the first and sometimes don't!!
                        if c == 1 && b == 1
                            % So firt block of each session, condition 11,
                            % are the ones that were modified manually
                            % (none others with new configuration). So if
                            % it's an 11_b1 (or whichever block is first in
                            % session_block_array), then continue to next
                            % block
                            continue
                        end
                        results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        infolder = find(results);
                        if isempty(infolder) % if no coincidence exists
                            continue 
                        end
                        if size(infolder,2)> 1 % in case of more than one coincidence (should not happen)
                           error(['more than one folder for ' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        end

                        % Get target channel file name
                        target_channel_name = [root_dir_bs '/data/' participant{p} '/' folders(infolder).name '/channel_vectorview306_acc1.mat'];

                        % Apply each of the transformations found in the original file
                        for lt = 1:length(x) % for each of the transformations found
                            % Select column of transformation
                            column_transf = y(lt); % Meg always
                            % Meg always, because EEG can have the 0 and 1 matrix when projecting electrodes to surface
                            Transf = TransfMeg{1,column_transf};
                            % Apply transformation (channel_apply_transf is a brainstorm function)
                            ChannelMats = channel_apply_transf(target_channel_name,  Transf);
                            % It's automatically saved
                        end
                    end
                end
            end      
            
            % Second, IN ALL blocks (including first), project EEG electrodes
            % to surface

            folders = dir([root_dir_bs '/data/' participant{p} '/']);
            for w = 1:length(wave) % always gonna be LLR now
                for c = 1:length(condition) % it's going to have to be applied to all 15 condition folders
                    for b = 1:length(block) % ALSO FIRST BLOCK NEEDS EEG PROJECTION TO SURFACE
                        results = contains({folders.name},[wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        infolder = find(results);
                        if isempty(infolder) % if no coincidence exists
                            continue 
                        end
                        if size(infolder,2)> 1 % in case of more than one coincidence (should not happen)
                           error(['more than one folder for ' wave{w} '_' condition{c} session{s} '_' block{b} '_normal']);
                        end

                        % Get name of any brainstorm file in the folder (from which
                        % channel file will be retrieved by function)
                        dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' folders(infolder).name]);
                        sweeps_norm = contains({dir_sweeps.name},'_trial');
                        sweeps_list = find(sweeps_norm);
                        position = sweeps_list(1);  % first trial found
                        file_name = [participant{p} '/' folders(infolder).name '/' dir_sweeps(position).name];
                        sFiles = file_name;

                        % Project electrodes to surface
                        sFiles = bst_process('CallProcess', 'process_channel_project', sFiles, []);
                    end
                end
            end
            
        end
    end
end

save([root_dir '/no_manual_transf.mat'],'no_manual_transf');

clearvars('-except', initialVars{:});
disp 'DONE WITH APLYING TRANSFORMATION MATRICES TO OTHER BLOCKS WITHIN SESSION (LLR)!!!'
disp(datetime)
toc
