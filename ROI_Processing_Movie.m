%% Actual fixation positions -> fixation information relative to face
close all
clear all
clc

%% load movie fixes & face box info
load('/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Analysis/Movie_Fixes')
movienames_run1 = {'MaisieWilliams', 'AmyPoehler', 'GordonRamsey', 'DanaCarvey', 'DonaldGlover', 'JohnMulaney', 'JohnOliver', 'NickOfferman', 'OprahWinfrey', 'Retta'};
movienames_run2 = {'WayneGretzky', 'Stanley', 'LucyLiu', 'TrevorNoah', 'AubreyPlaza', 'DaveedDiggs', 'KristenStewart', 'BobbyFlay', 'AdamScott', 'WyattCenac'};
movienames_run3 = {'MargaretBrennan', 'JonStewart', 'JJAbrams', 'MaxGreenfield', 'Ellen', 'BenedictCumberbatch', 'TinaFey', 'EllieKemper', 'AlexisOhanian', 'EmilyBlunt'}; 
movienames_run4 = {'MilaKunis', 'YaraShahidi', 'LillySingh', 'KristenBell', 'JohnKrasinski', 'AishaTaylor', 'KumailNanjiani', 'RyanReynolds', 'BarackObama', 'AngelaDuckworth'};

movienames_all = {movienames_run1; movienames_run2; movienames_run3; movienames_run4};

fix_mat = table2array(fix_report_tot);
n_subs = numel(unique(fix_mat(:,1)));
fps = 29.97; % frames per second
spf = 1/fps; % seconds per frame
mf = 0.10; % mouth line factor

%%
% loop through each subject
% loop through each run
% loop through each movie
% --> load in face info for each from of selected movie
% loop through each fix in that movie
% --> get relevant frame for fix
% --> determine if fix is on face (=0), upper face (=1), lower face(=2) 

% Prepare final var matrices
region_by_fix = [];
movie_avgd = [];
run_avgd = [];
sub_avgd = [];

for i = 1:n_subs
    % get data for current sub
    sub_data = fix_mat(fix_mat(:,1) == i,:);
    n_runs = numel(unique(sub_data(:,2)));
    
    for j = 1:n_runs
        % get data for current run
        run_data = sub_data(sub_data(:,2) == j,:);
        n_movies = numel(unique(run_data(:,3)));
        
        for m = 1:n_movies
            movie_num = m;
                                    
            % get data for current movie
            movie_data = run_data(run_data(:,3) == m,:); % fix info for this subject for this movie
            c_movie_str = char(movienames_all{j}(m)); % name of current movie
            n_fixes = numel(movie_data(:,1)); 
            
            % load in face box information from detection script (x y width height) 
            load(['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Analysis/Adj_FacePos/' c_movie_str '_adj.mat'])
                        
            %
            mean_size = round(nanmean(box_pos(:,3)));
            n_frame = length(box_pos);
            total_t = (n_frame-1)*spf; % in seconds
            timings = (0:spf:total_t)'; % onsets of frames in ms
            frame_at_time = [timings box_pos];
            
            for k = 1:n_fixes                    
                % get location of fixation
                fix_data = movie_data(k,:);
                fix_x = movie_data(k,4);
                fix_y = movie_data(k,5);
                
                % get timing of fixation (convert ms -> seconds)
                fix_start = movie_data(k,7)/1000;
                fix_end = movie_data(k,8)/1000;
                
                % which frame(s) were seen during this fix?
                frames_after_start = frame_at_time(frame_at_time(:,1) > fix_start, :);
                frames_during = frames_after_start(frames_after_start(:,1) < fix_end, :);
                n_frames_seen_during_c_fix = length(frames_during(:,1)); % number of frames
                
                % restarts region counter for this fix
                region = [];
                
                % runs through all frames shown during this fix
                for f = 1:n_frames_seen_during_c_fix 

                    c_face = frames_during(f,2:5);
                    
                    % where is the face in this frame? (x, y pos)
                    frame_coords = [c_face(1) c_face(2); % top_left
                        c_face(1)+c_face(3) c_face(2); % top right
                        c_face(1) c_face(2)+c_face(4); % bottom left
                        c_face(1)+c_face(3) c_face(2)+c_face(4)]; % bottom right
                    
                    frame_coords = frame_coords * 1.5;
                    half_dist = (frame_coords(3,2) - frame_coords(1,2))/2;
                    halfway_split = frame_coords(1,2) + half_dist + (mf*mean_size);

                    % is the fix on the face in this frame?
                    if fix_x > frame_coords(1,1) && fix_x < frame_coords(2,1) && fix_y > frame_coords(1,2) && fix_y < frame_coords(3,2)
                        if fix_y > halfway_split % lower half = 1; used > because (0,0) is at top left
                            region(f) = 1;
                        elseif fix_y < halfway_split % upper half = 2; used < because (0,0) is at top left
                            region(f) = 2;
                        end
                    else 
                        region(f) = 0; % not on face
                    end
                end
                
                % Calculate % in each region
                total_frames = length(region);
                p_lower = sum(region==1)/total_frames;
                p_upper = sum(region==2)/total_frames;
                p_notface = sum(region==0)/total_frames;
                
                % line in fix_report_tot: position + % time by region + frames seen
                region_by_fix = [region_by_fix; fix_data p_lower p_upper p_notface n_frames_seen_during_c_fix];
                
            end
            
            % summarize over sub, run, movie to get avg for each movie
            c_sub_data = region_by_fix(region_by_fix(:,1) == i,:);
            c_run_data = c_sub_data(c_sub_data(:,2) == j,:);
            c_movie_data = c_run_data(c_run_data(:,3) == m,:);
            
            
            if c_movie_data(:,9:12) ==  
            else
            movie_avgd = [movie_avgd; i j m wmean(c_movie_data(:,4:5), c_movie_data(:,12)) wmean(c_movie_data(:,9:11), c_movie_data(:,12))];
            end
        end
        
        run_avgd = [run_avgd; i j wmean(c_run_data(:,4:5), c_run_data(:,12)), wmean(c_run_data(:,9:11), c_run_data(:,12))];
        
    end
    
    sub_avgd = [sub_avgd; i wmean(c_sub_data(:,4:5), c_sub_data(:,12)), wmean(c_sub_data(:,9:11), c_sub_data(:,12))];
    
end

%% Save
display('Save...');

head_line = {'Sub' 'Run' 'Movie' 'x' 'y' 'FixDur' 'FixStart' 'FixEnd', 'P_Lower', 'P_Upper', 'P_NotFace', 'NFrames'};
region_report = array2table(region_by_fix,'VariableNames',head_line);

head_line = {'Sub' 'Run' 'Movie' 'Avg_x' 'Avg_y','P_Lower', 'P_Upper', 'P_NotFace'};
movie_report_tot = array2table(movie_avgd,'VariableNames',head_line);

head_line = {'Sub' 'Run' 'Avg_x' 'Avg_y','P_Lower', 'P_Upper', 'P_NotFace'};
run_report_tot = array2table(run_avgd,'VariableNames',head_line);

head_line = {'Sub' 'Avg_x' 'Avg_y','P_Lower', 'P_Upper', 'P_NotFace'};
sub_report_tot = array2table(sub_avgd,'VariableNames',head_line);

save('MovieFixes_ROI_Tab.mat','fix_report_tot', 'movie_report_tot', 'run_report_tot', 'sub_avgd');

display('Analysis done!');
