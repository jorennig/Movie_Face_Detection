%% Takes eye movement data, makes fixation matrix with one line for each frame
% ascii file (converted from edf) --> 
% matrix: [sub, movie, frame, x, y, fix duration]

%% PREP
clear all
close all
clc

if exist('Movie_Fixes.mat','file') == 2
    
    check_data = 1;
    
    load('Movie_Fixes.mat');
    
    nb_sub_data = unique(table2array(fix_report_tot(:,1))); % Number subjects already processed
    
    % Get all asc files
    result_files = dir('/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_ET/*Movie*.asc');
    result_files = {result_files(:).name}';
    
    % Get number of subjects not yet processed
    sub_files = zeros(length(result_files),1);
    
    for i = 1:length(result_files)
        
        sub_c = result_files{i};
        sub_files(i) = str2double(sub_c(1:2));
        
    end
    
    nb_sub_files = unique(sub_files);
    
    % Difference between data and files
    diff_files_data = setdiff(nb_sub_files,nb_sub_data);
    
    if isempty(diff_files_data) % in case all raw eye tracking data is already processed terminate script
        display('Nothing to do!');
        return
    end

    idx = find(sub_files >= min(diff_files_data));
    result_files = result_files(idx);
    fix_report_tot_old = table2array(fix_report_tot);
else
    
    check_data = 0;
    
    % Get all asc files
    result_files = dir('/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_ET/*Movie*.asc');
    result_files = {result_files(:).name}';
    
end

% Prepare data matrices
fix_report = [];
fix_report_tot = [];

%% LOAD DATA
for i = 1:numel(result_files)
    
    full_file_path = ['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_ET/' result_files{i}];
    data_samples_run = zeros(1,5);

    ET = strsplit(fileread(full_file_path), '\n')';
    
    result_file_name = result_files{i};
    
    sub_num = result_file_name(1:2);
    sub_num = str2double(sub_num);
    
    % Get run
    idxe = strfind(result_file_name,'_');
    result_run = result_file_name(idxe(1)+1:idxe(2)-1);
    run_num = str2double(result_run(4)); % Number of runs
    
    % Get sampling rate
    idx_sr = strfind(ET,'RECCFG');
    idx_sr = find(not(cellfun('isempty', idx_sr)));
    idx_sr = ET(idx_sr(1));
    idx_sr = cellfun(@(x) strsplit(x,' '),idx_sr(:),'uni',0);
    idx_sr = idx_sr{1,1};
    sampling_rate = str2double(cell2mat(idx_sr(4)));
    
    % Get start of experiment
    idx = strfind(ET,'SYNCTIME');
    idx_sync = find(not(cellfun('isempty', idx)));
    ET = ET(idx_sync:end);
    
    idx = strfind(ET,'start movie'); % Movei start & end
    idx_mov_start = find(not(cellfun('isempty', idx)));
    idx = strfind(ET,'end movie');
    idx_mov_end = find(not(cellfun('isempty', idx)));    
    
    fprintf('-- Subject %d, RUN %d --\n',sub_num,run_num)

    trial_report = []; % Pre-define trial summary
    fix_report = []; % Pre-define fixation summary

    %% Trialwise processing movies
    for j = 1:length(idx_mov_start) % Movie data
        
        fprintf('Movie Nr %d of %d\n',j,length(idx_mov_start))
        data_trial = ET(idx_mov_start(j)+1:idx_mov_end(j)-1); % All data from trial
        
        % Get movie ID
        movie_idn = [];
        idx_t = 0;
        
        while isempty(movie_idn)
            movie_idv = ET(idx_mov_start(j)-idx_t); % Read out MOVIEID
            movie_idc = strfind(movie_idv,'MOVIEID');
            movie_idn = find(not(cellfun('isempty', movie_idc)));
            idx_t = idx_t + 1;
        end

        movie_idl = movie_idv(movie_idn);
        movie_idstr = movie_idl{1};
        movie_idp = strfind(movie_idl,'MOVIEID');
        movie_idp = movie_idp{1};
        movie_id = movie_idstr(movie_idp+8:end);
        movie_id = str2double(movie_id); % Movie ID/number        
        
        movie_start_time = str2double(movie_idstr(movie_idp-8:movie_idp-2))+1;
        
        name_trialc = ET(idx_mov_start(j)); % Read out movie name
        name_trialstr = name_trialc{1};                
        name_trials = strfind(name_trialstr,'movie');
        name_triale = strfind(name_trialstr,'.');
        name_trial = name_trialstr(name_trials+8:name_triale-1); % Movie name        
        
        %% Process data
        data_trial(1:100) = []; % Remove first few samples since they are part of saccade prior to first fixation
        
        num_sacc_trial_idx = strfind(data_trial,'SSACC'); % Number of saccades
        num_sacc_trial = numel(find(not(cellfun('isempty', num_sacc_trial_idx))));
        sacc_trial = find(not(cellfun('isempty', num_sacc_trial_idx))); % Start saccade
        num_sacc_trial = num_sacc_trial - 1; % remove first saccade after trial onset
        
        num_blink_trial_idx = strfind(data_trial,'SBLINK'); % Number of blinks
        blink_trial = find(not(cellfun('isempty', num_blink_trial_idx))); % Stat blink
        num_blink_trial = numel(find(not(cellfun('isempty', num_blink_trial_idx))));            
        
        % Check if we have Start and End marks for fixations
        % Get all fixation starts
        fix_start_idxv = strfind(data_trial,'SFIX'); % Starts of fixations
        fix_start_idx = find(not(cellfun('isempty', fix_start_idxv)));

        % Get all fixation ends
        fix_end_idxv = strfind(data_trial,'EFIX'); % Ends of fixations
        fix_end_idx = find(not(cellfun('isempty', fix_end_idxv)));

%         % In case fixation marks lie outside movie marks
%         if isempty(fix_start_idx) == 1 && numel(fix_end_idx) == 1
%             fix_start_idx = 2;
%         end
%         if isempty(fix_end_idx) == 1 && numel(fix_start_idx) == 1
%             fix_end_idx = length(data_trial)-1;
%         end
%         if isempty(fix_end_idx) == 1 && isempty(fix_start_idx) == 1
%             fix_start_idx = 2;
%             fix_end_idx = length(data_trial)-1;
%         end
% 
%         % In case last fixation start is at the very end of data vector
%         if fix_start_idx(end) > length(data_trial)-20
%             fix_start_idx(end) = [];            
%         end
% 
%         % In case first fixation end is at the very beginning of data vector
%         if fix_end_idx(1) < 20
%             fix_start_idx(end) = [];
%         end

        % If we can find fixation start and end points
        sample_thresh = 200;
        
        if isempty(fix_start_idx) == 0 && isempty(fix_end_idx) == 0
            
            % If we do NOT have an identical numbers of fixation starts and ends
            if numel(fix_start_idx) ~= numel(fix_end_idx)
                
                % If we have more fixation starts than ends
                % We only have an additional start at the end of the trial
                if numel(fix_start_idx) > numel(fix_end_idx)
                    
                    if fix_start_idx(end) < length(data_trial)-sample_thresh
                        fix_end_idx = [fix_end_idx;length(data_trial)];
                    elseif fix_start_idx(end) >= length(data_trial)-sample_thresh
                        fix_start_idx(end) = [];
                    end
                    
                % If we have more fixation ends than starts
                % We only have an additional end at the start of the trial
                elseif numel(fix_start_idx) < numel(fix_end_idx)
                    
                    if fix_end_idx(1) > sample_thresh
                        fix_start_idx = [2; fix_start_idx];
                    elseif fix_end_idx(1) <= sample_thresh
                        fix_end_idx(1) = [];
                    end
                    
                end
            
            % If we DO have an identical number of fixation starts and ends
            % but first fixation END occurs earlier than first START
            elseif numel(fix_start_idx) == numel(fix_end_idx)
                
                % If fixations are cut off at both ends
                if fix_end_idx(1) < fix_start_idx(1)
                    
                    comp = [fix_start_idx fix_end_idx];
                    
                    % If fixations at both ends are long enough to be used
                    if fix_start_idx(1) > sample_thresh && fix_end_idx(end) < length(data_trial)-sample_thresh
                        fix_start_idx = [2; fix_start_idx];
                        fix_end_idx = [fix_end_idx; length(data_trial)-1];
                    
                    % If fixations at both ends are too short to be used
                    elseif fix_start_idx(1) <= sample_thresh && fix_end_idx(end) >= length(data_trial)-sample_thresh
                        fix_start_idx(end) = [];
                        fix_end_idx(1) = [];
                        
                    % If fixation at the start is long enough, too short at
                    % the end
                    elseif fix_start_idx(1) > sample_thresh && fix_end_idx(end) >= length(data_trial)-sample_thresh
                        fix_start_idx(end) = [];
                        fix_start_idx = [2; fix_start_idx];
                    
                    % If fixation at the end is long enough, too short at
                    % the start
                    elseif fix_start_idx(1) <= sample_thresh && fix_end_idx(end) < length(data_trial)-sample_thresh
                        fix_end_idx(1) = [];
                        fix_end_idx = [fix_end_idx; length(data_trial)-1];
                    end
                    
                end
                
            end
            
        else
            fix_start_idx = 2;
            fix_end_idx = length(data_trial) - 2;
        end
                                
        % Pre-allocate for all fixations in trial
        et_data_list = [];

        %% Process all samples that are fixations
        for p = 1:numel(fix_start_idx)
            
            data_fix = data_trial(fix_start_idx(p)+1:fix_end_idx(p)-1);           
            data_fixv = cellfun(@(x) strsplit(x,' '),data_fix(:),'uni',0);
            
            rel = 2; % relevant values from data_fix
            data_fixn = zeros(numel(data_fixv),rel);
            
            start = data_fixv{1};
            start_time = str2double(cell2mat(start(1))) - movie_start_time;
            
            for q = 1:numel(data_fixv)
                dat_cur = data_fixv{q};
                d_time(q,:) = str2double(cell2mat(dat_cur(1)));
                d_x = str2double(cell2mat(dat_cur(2)));
                d_y = str2double(cell2mat(dat_cur(3)));
                data_fixn(q,:) = [d_x d_y];
            end
            
            end_time = d_time(numel(data_fixv)) - movie_start_time;
            
            if size(data_fixn,1) == 1
                fix_mean = data_fixn;
            else
                fix_mean = nanmean(data_fixn);
            end
            
            fix_duration = size(data_fixn,1)/sampling_rate*1000;
            
            % Fixation means and durations in list
            et_data_list(p,:) = [sub_num run_num movie_id fix_mean end_time-start_time start_time end_time];
            
        end
        
        % Save means to fixation report
        fix_report = [fix_report; et_data_list];
        
    end

    % Save to overall fixation report        
    fix_report_tot = [fix_report_tot; fix_report];
    
end

%% SAVE
head_line = {'Sub' 'Run' 'Movie' 'x' 'y' 'FixDur' 'FixStart' 'FixEnd'};

display('Save...');

if check_data == 1 % in case of add-on analysis
    fix_report_tot = [fix_report_tot_old; fix_report_tot];
    fix_report_tot = array2table(fix_report_tot,'VariableNames',head_line);
    %trial_report_tot = [trial_report_tot_old; trial_report_tot];
    save('Movie_Fixes.mat','fix_report_tot');
else
    fix_report_tot = array2table(fix_report_tot,'VariableNames',head_line);
    save('Movie_Fixes.mat','fix_report_tot');
end

display('Analysis done!');
