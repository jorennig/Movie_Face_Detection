%% Clean up
close all
clear all
clc

%% Plot variables
lw = 1.1; % Line Width
bw = 0.8; % Barwidth
fsb = 16; % Font Size big
fsbb = 18; % Font Size small
path = [pwd '/Heat_Maps']; % Current path
rez = 300; % resolution (dpi) of final graphic
format = '-djpeg'; % image format

area = 20;
mf = 0.10; % mouth line factor

%% Load data and get face box sizes
load('/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Analysis/Movie_Fixes.mat')
fix_mat = table2array(fix_report_tot);

% Correct run1 for subject 1
run1 = fix_mat(fix_mat(:,2)==1,:);
run1(:,7:8) = run1(:,7:8) - 10000000;
fix_mat(fix_mat(:,2)==1,:) = run1;

n_subs = numel(unique(fix_mat(:,1)));
fps = 29.97; % frames per second
spf = 1/fps; % seconds per frame
monitor_res = [1920,1080];
display_center = monitor_res/2;
img_res = [1280,720];
img_center = img_res/2;

load('/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Analysis/Movie_Fixes.mat')
movienames_run1 = {'MaisieWilliams', 'AmyPoehler', 'GordonRamsey', 'DanaCarvey', 'DonaldGlover', 'JohnMulaney', 'JohnOliver', 'NickOfferman', 'OprahWinfrey', 'Retta'};
movienames_run2 = {'WayneGretzky', 'Stanley', 'LucyLiu', 'TrevorNoah', 'AubreyPlaza', 'DaveedDiggs', 'KristenStewart', 'BobbyFlay', 'AdamScott', 'WyattCenac'};
movienames_run3 = {'MargaretBrennan', 'JonStewart', 'JJAbrams', 'MaxGreenfield', 'Ellen', 'BenedictCumberbatch', 'TinaFey', 'EllieKemper', 'AlexisOhanian', 'EmilyBlunt'}; 
movienames_run4 = {'MilaKunis', 'YaraShahidi', 'LillySingh', 'KristenBell', 'JohnKrasinski', 'AishaTaylor', 'KumailNanjiani', 'RyanReynolds', 'BarackObama', 'AngelaDuckworth'};

movienames_all = {movienames_run1, movienames_run2; movienames_run3, movienames_run4};
movienames_list = [movienames_run1, movienames_run2; movienames_run3, movienames_run4];

% Calculating average box size to use
box_size = 400;

for i = 1:numel(movienames_list)

    c_movie_str = char(movienames_list{i});

    % load in face box information from detection script (x y width
    % height) & adjust to monitor size
    load(['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Analysis/Adj_FacePos/' c_movie_str '_adj.mat'])
    
    box_size_tot(i,:) = box_pos(1,3:4);
    box_pos_tot(i,:) = round(nanmean(box_pos(:,1:2)));
    
end

mean_box_size = round(mean(box_size_tot));

factor = box_size_tot(:,1)/box_size;
factor_run = reshape(factor,[10,4]);

%% Visualize original data per subject/run (fixations)
sub = unique(fix_mat(:,1));

for i = 1:numel(sub)
    
    fprintf('-- Subject %d --\n',sub(i))

    data_s = fix_mat(fix_mat(:,1)==sub(i),:);
    fix_overall = sum(data_s(:,3));
    
    data_s(any(isnan(data_s),2),:) = [];
    runs = unique(data_s(:,2));
    
    plot_idx = 0;
    box_count = 0;

    for j = 1:numel(runs)
        
        data_r = data_s(data_s(:,2)==runs(j),:);
        movies = unique(data_r(:,3));
        
        for l = 1:numel(movies)
            
            data_m = data_r(data_r(:,3)==movies(l),:);
            
            % Plot
            plot_idx = plot_idx + 1;
%             if plot_idx == 8
%                 plot_idx = plot_idx + 1;
%             end
            subplot(numel(runs),numel(movies),plot_idx)
            scatter(data_m(:,4),data_m(:,5))
            
            hold on
            box_count = box_count + 1;
            box_size_c = box_size_tot(box_count);
            box_pos_c = box_pos_tot(box_count,:);
            
            box_coords = [box_pos_c(1) box_pos_c(2); % top_left
                box_pos_c(1)+box_size_c box_pos_c(2); % top right
                box_pos_c(1) box_pos_c(2)+box_size_c; % bottom left
                box_pos_c(1)+box_size_c box_pos_c(2)+box_size_c]; % bottom right

            box_coords = round(box_coords * 1.5);
            box_size_c_n = box_coords(3,2)-box_coords(2,2);
            split_line = round(box_coords(1,2) + (box_size_c_n/2) + (mf*box_size_c_n));
            
            line([box_coords(1,1),box_coords(2,1)],[box_coords(1,2) box_coords(2,2)],'LineWidth',2,'Color','k'); % Line horizontal low
            line([box_coords(1,1),box_coords(2,1)],[box_coords(3,2) box_coords(4,2)],'LineWidth',2,'Color','k'); % Line horizontal high
            line([box_coords(1,1),box_coords(1,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical left
            line([box_coords(2,1),box_coords(2,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical right
            line([box_coords(1,1),box_coords(2,1)],[split_line split_line],'LineWidth',2,'Color','k'); % Split line
            
            set(gca,'Ydir','reverse');
            
            xlim([400 1600]);
            ylim([0 1000]);
            
            title(['Movie ' num2str(l)]);
            
        end
    end
    
    % Save plot
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'Units','centimeters','Position',[10 20 80 15]);
    f = gcf; % f is the handle of the figure you want to export
    if sub(i) < 10
        name = ['Fixations_org_Mov_Sub0' num2str(sub(i))];
    else
        name = ['Fixations_org_Mov_Sub' num2str(sub(i))];
    end

    set(gcf,'PaperPositionMode','auto')
    print(f,fullfile(path,name),format,['-r',num2str(rez)])

    close all

end

%% Transform data
data_corr = [];
for i = 1:n_subs
    % pick data for current sub
    sub_data = fix_mat(fix_mat(:,1) == i,:);
    n_runs = numel(unique(sub_data(:,2)));
    
    sub_data_corr = [];
    for j = 1:n_runs
        % pick data for current run
        run_data = sub_data(sub_data(:,2) == j,:);
        n_movies = numel(unique(run_data(:,3)));
        
        run_data_corr = [];
        for m = 1:n_movies
            movie_num = m;
                                    
            % get data for current movie
            movie_data = run_data(run_data(:,3) == m,:); % fix info for this subject for this movie
            c_movie_str = char(movienames_all{j}(m));
            n_fixes = numel(movie_data(:,1));
            
            % load in face box information from detection script (x y width
            % height) & adjust to monitor size
            load(['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Analysis/Adj_FacePos/' c_movie_str '_adj.mat'])
                        
            %
            mean_size = round(nanmean(box_pos(:,3)));
            n_frame = length(box_pos);
            total_t = (n_frame-1)*spf; % in seconds
            timings = (0:spf:total_t)'; % onsets of frames in ms
            frame_at_time = [timings box_pos];
            avg_position_frames_vid = mean(frame_at_time(:,2:3));
            avg_position_frames_vid_center = avg_position_frames_vid + (mean_size/2);
            diff_avg_frames_center_disp_center = img_center - avg_position_frames_vid_center;
            
            diff_position_frames_vid_fix = [];
            for k = 1:n_fixes
                
                % get location of fixation
                fix_data = movie_data(k,:);
                fix_x = movie_data(k,4);
                fix_y = movie_data(k,5);
                
                % get timing of fixation (in seconds)
                fix_start = movie_data(k,7)/1000;
                fix_end = movie_data(k,8)/1000;
                
                % which frame(s) were seen during this fix?
                frames_after_start = frame_at_time(frame_at_time(:,1) > fix_start, :);
                frames_during = frames_after_start(frames_after_start(:,1) < fix_end, :);
                %n_frames_seen_during_c_fix = length(frames_during(:,1)); % number of frames
                
                avg_position_frames_fix = mean(frames_during(:,2:3));
                diff_position_frames_vid_fix(k,:) = avg_position_frames_vid - avg_position_frames_fix;
                                
            end
            
            % Subtract differences of mean reference frame from frames
            % associated with respective fixation
            movie_data_corr = movie_data(:,4:5) + diff_position_frames_vid_fix;
            
            % Center the data to center of the monitor
            movie_data_corr = movie_data_corr + diff_avg_frames_center_disp_center;
            
            % Scale coordinates according to box size
            fac_c = factor_run(m,j);
            %movie_data_corr = round(movie_data_corr*fac_c);
            movie_data_corr = round(movie_data_corr*1);
            
            movie_data_corr = [movie_data(:,1:3), movie_data_corr, movie_data(:,6)];
            
            run_data_corr = [run_data_corr;movie_data_corr];
            
        end
        
        sub_data_corr = [sub_data_corr;run_data_corr];
        
    end
    
    data_corr = [data_corr;sub_data_corr];
    
end

%% Visualize transformed data per subject (fixations)
close all

for i = 1:numel(sub)
    
    fprintf('-- Subject %d --\n',sub(i))

    data_s = data_corr(data_corr(:,1)==sub(i),:);
    fix_overall = sum(data_s(:,3));
    
    data_s(any(isnan(data_s),2),:) = [];
    runs = unique(data_s(:,2));
    
    plot_idx = 0;
    for j = 1:numel(runs)
        
        data_r = data_s(data_s(:,2)==runs(j),:);
        movies = unique(data_r(:,3));
        
        for l = 1:numel(movies)
            
            data_m = data_r(data_r(:,3)==movies(l),:);
            
            % Plot
            plot_idx = plot_idx + 1;

            subplot(numel(runs), numel(movies),plot_idx)
            scatter(data_m(:,4),data_m(:,5))
            
            hold on
            
            box_pos_c = [img_center(1)-box_size/2 img_center(2)-box_size/2];
            box_coords = [box_pos_c(1) box_pos_c(2); % top_left
                box_pos_c(1)+box_size box_pos_c(2); % top right
                box_pos_c(1) box_pos_c(2)+box_size; % bottom left
                box_pos_c(1)+box_size box_pos_c(2)+box_size]; % bottom right

            box_coords = round(box_coords * 1.5);
            box_size_c_n = box_coords(3,2)-box_coords(2,2);
            split_line = round(box_coords(1,2) + (box_size_c_n/2) + (mf*box_size_c_n));
            
            line([box_coords(1,1),box_coords(2,1)],[box_coords(1,2) box_coords(2,2)],'LineWidth',2,'Color','k'); % Line horizontal low
            line([box_coords(1,1),box_coords(2,1)],[box_coords(3,2) box_coords(4,2)],'LineWidth',2,'Color','k'); % Line horizontal high
            line([box_coords(1,1),box_coords(1,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical left
            line([box_coords(2,1),box_coords(2,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical right
            line([box_coords(1,1),box_coords(2,1)],[split_line split_line],'LineWidth',2,'Color','k'); % Split line
            
            set(gca,'Ydir','reverse');

            xlim([400 1400]);
            ylim([0 1200]);
                        
            title(['Movie ' num2str(l)]);

        end
    end
    
    % Save plot
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'Units','centimeters','Position',[10 20 80 15]);
    f = gcf; % f is the handle of the figure you want to export
    if sub(i) < 10
        name = ['Fixations_trans_Mov_Sub0' num2str(sub(i))];
    else
        name = ['Fixations_trans_Mov_Sub' num2str(sub(i))];
    end

    set(gcf,'PaperPositionMode','auto')
    print(f,fullfile(path,name),format,['-r',num2str(rez)])

    close all

end

%% Visualize transformed data per subject/run (heat maps)
sub = unique(data_corr(:,1));

for i = 1:numel(sub)
    
    fprintf('-- Subject %d --\n',sub(i))

    data_s = data_corr(data_corr(:,1)==sub(i),:);
    fix_overall = sum(data_s(:,3));
    
    data_s(data_s(:,4) > monitor_res(1)-area,:) = NaN;
    data_s(data_s(:,5) > monitor_res(2)-area,:) = NaN;

    data_s(data_s(:,4) < area,:) = NaN;
    data_s(data_s(:,5) < area,:) = NaN;

    data_s(any(isnan(data_s),2),:) = [];
    runs = unique(data_s(:,2));
    
    plot_idx = 0;
    for j = 1:numel(runs)
        
        data_r = data_s(data_s(:,2)==runs(j),:);
        movies = unique(data_r(:,3));
        
        for l = 1:numel(movies)
            
            data_m = data_r(data_r(:,3)==movies(l),:);
            
            % Heat map template
            heat_map = zeros(monitor_res(2),monitor_res(1));

            % Square area
            for k = 1:size(data_m,1)
                
                val = data_m(k,6);
                heat_map(data_m(k,5)-area:data_m(k,5)+area,data_m(k,4)-area:data_m(k,4)+area) = heat_map(data_m(k,5)-area:data_m(k,5)+area,data_m(k,4)-area:data_m(k,4)+area) + val;

            end

            % Calculate percent
            heat_map = heat_map/fix_overall*100;

            % Smooth heat map
            g_kernel = fspecial('gaussian',60,60);
            heat_map = filter2(g_kernel,heat_map);

            % Set threshold
            hm_values = sort(reshape(heat_map,[1 size(heat_map,1)*size(heat_map,2)]));
            hm_values(hm_values==0) = [];

            thr_pc = 70; % percent of lower values to be cut off
            vals = unique(hm_values);
            vals_sum = numel(vals);
            val_crit = round(vals_sum*(thr_pc/100));
            cut_off = vals(val_crit);

            % Apply threshold
            heat_map(heat_map < cut_off) = 0;

            %heat_map = heat_map(monitor_res(2)/2-stim_size(2)/2:monitor_res(2)/2+stim_size(2)/2-1,monitor_res(1)/2-stim_size(1)/2:monitor_res(1)/2+stim_size(1)/2-1);

            % Plot
            plot_idx = plot_idx + 1;
%             if plot_idx == 8
%                 plot_idx = plot_idx + 1;
%             end
            subplot(4,10,plot_idx)
            cm_c = jet;
            cm_c(1,:) = [1 1 1];
            colormap(cm_c);
            imagesc(heat_map);

            hold on
            box_pos_c = [img_center(1)-box_size/2 img_center(2)-box_size/2];
            box_coords = [box_pos_c(1) box_pos_c(2); % top_left
                box_pos_c(1)+box_size box_pos_c(2); % top right
                box_pos_c(1) box_pos_c(2)+box_size; % bottom left
                box_pos_c(1)+box_size box_pos_c(2)+box_size]; % bottom right

            box_coords = round(box_coords * 1.5);
            box_size_c_n = box_coords(3,2)-box_coords(2,2);
            split_line = round(box_coords(1,2) + (box_size_c_n/2) + (mf*box_size_c_n));
            
            line([box_coords(1,1),box_coords(2,1)],[box_coords(1,2) box_coords(2,2)],'LineWidth',2,'Color','k'); % Line horizontal low
            line([box_coords(1,1),box_coords(2,1)],[box_coords(3,2) box_coords(4,2)],'LineWidth',2,'Color','k'); % Line horizontal high
            line([box_coords(1,1),box_coords(1,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical left
            line([box_coords(2,1),box_coords(2,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical right
            line([box_coords(1,1),box_coords(2,1)],[split_line split_line],'LineWidth',2,'Color','k'); % Split line
            
            %set(gca,'xtick',[]);
            %set(gca,'ytick',[]);
            %box off
            xlim([400 1400]);
            ylim([0 1200]);
            
            title(['Movie ' num2str(l)]);

        end
    end
    
    % Save plot
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'Units','centimeters','Position',[10 20 80 15]);
    f = gcf; % f is the handle of the figure you want to export
    if sub(i) < 10
        name = ['Heatmap_trans_Mov_Sub0' num2str(sub(i))];
    else
        name = ['Heatmap_trans_Mov_Sub' num2str(sub(i))];
    end

    set(gcf,'PaperPositionMode','auto')
    print(f,fullfile(path,name),format,['-r',num2str(rez)])

    close all

end

%% Visualize transformed data per run (fixations)
close all

for i = 1:numel(sub)
    
    fprintf('-- Subject %d --\n',sub(i))

    data_s = data_corr(data_corr(:,1)==sub(i),:);
    fix_overall = sum(data_s(:,3));
    
    data_s(any(isnan(data_s),2),:) = [];
    runs = unique(data_s(:,2));
    
    for j = 1:numel(runs)
        
        data_r = data_s(data_s(:,2)==runs(j),:);
        
        subplot(1,4,j)
        scatter(data_r(:,4),data_r(:,5))

        hold on

        box_pos_c = [img_center(1)-box_size/2 img_center(2)-box_size/2];
        box_coords = [box_pos_c(1) box_pos_c(2); % top_left
            box_pos_c(1)+box_size box_pos_c(2); % top right
            box_pos_c(1) box_pos_c(2)+box_size; % bottom left
            box_pos_c(1)+box_size box_pos_c(2)+box_size]; % bottom right

        box_coords = round(box_coords * 1.5);
        box_size_c_n = box_coords(3,2)-box_coords(2,2);
        split_line = round(box_coords(1,2) + (box_size_c_n/2) + (mf*box_size_c_n));

        line([box_coords(1,1),box_coords(2,1)],[box_coords(1,2) box_coords(2,2)],'LineWidth',2,'Color','k'); % Line horizontal low
        line([box_coords(1,1),box_coords(2,1)],[box_coords(3,2) box_coords(4,2)],'LineWidth',2,'Color','k'); % Line horizontal high
        line([box_coords(1,1),box_coords(1,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical left
        line([box_coords(2,1),box_coords(2,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical right
        line([box_coords(1,1),box_coords(2,1)],[split_line split_line],'LineWidth',2,'Color','k'); % Split line

        set(gca,'Ydir','reverse');

        xlim([400 1400]);
        ylim([0 1200]);

        title(['Run ' num2str(j)]);

    end
    
    % Save plot
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'Units','centimeters','Position',[10 20 40 15]);
    f = gcf; % f is the handle of the figure you want to export
    if sub(i) < 10
        name = ['Fixations_trans_Run_Sub0' num2str(sub(i))];
    else
        name = ['Fixations_trans_Run_Sub' num2str(sub(i))];
    end

    set(gcf,'PaperPositionMode','auto')
    print(f,fullfile(path,name),format,['-r',num2str(rez)])

    close all

end

%% Visualize transformed data per run (heat maps)
close all

for i = 1:numel(sub)
    
    fprintf('-- Subject %d --\n',sub(i))

    data_s = data_corr(data_corr(:,1)==sub(i),:);
    fix_overall = sum(data_s(:,3));
    
    data_s(data_s(:,4) > monitor_res(1)-area,:) = NaN;
    data_s(data_s(:,5) > monitor_res(2)-area,:) = NaN;

    data_s(data_s(:,4) < area,:) = NaN;
    data_s(data_s(:,5) < area,:) = NaN;

    data_s(any(isnan(data_s),2),:) = [];
    runs = unique(data_s(:,2));
    
    for j = 1:numel(runs)
        
        data_r = data_s(data_s(:,2)==runs(j),:);
        
        % Heat map template
        heat_map = zeros(monitor_res(2),monitor_res(1));

        % Square area
        for k = 1:size(data_r,1)

            val = data_r(k,6);
            heat_map(data_r(k,5)-area:data_r(k,5)+area,data_r(k,4)-area:data_r(k,4)+area) = heat_map(data_r(k,5)-area:data_r(k,5)+area,data_r(k,4)-area:data_r(k,4)+area) + val;

        end

        % Calculate percent
        heat_map = heat_map/fix_overall*100;

        % Smooth heat map
        g_kernel = fspecial('gaussian',60,60);
        heat_map = filter2(g_kernel,heat_map);

        % Set threshold
        hm_values = sort(reshape(heat_map,[1 size(heat_map,1)*size(heat_map,2)]));
        hm_values(hm_values==0) = [];

        thr_pc = 80; % percent of lower values to be cut off
        vals = unique(hm_values);
        vals_sum = numel(vals);
        val_crit = round(vals_sum*(thr_pc/100));
        cut_off = vals(val_crit);

        % Apply threshold
        heat_map(heat_map < cut_off) = 0;

        %heat_map = heat_map(monitor_res(2)/2-stim_size(2)/2:monitor_res(2)/2+stim_size(2)/2-1,monitor_res(1)/2-stim_size(1)/2:monitor_res(1)/2+stim_size(1)/2-1);

        % Plot
        subplot(1,4,j)
        cm_c = jet;
        cm_c(1,:) = [1 1 1];
        colormap(cm_c);
        imagesc(heat_map);

        hold on
        box_pos_c = [img_center(1)-box_size/2 img_center(2)-box_size/2];
        box_coords = [box_pos_c(1) box_pos_c(2); % top_left
            box_pos_c(1)+box_size box_pos_c(2); % top right
            box_pos_c(1) box_pos_c(2)+box_size; % bottom left
            box_pos_c(1)+box_size box_pos_c(2)+box_size]; % bottom right

        box_coords = round(box_coords * 1.5);
        box_size_c_n = box_coords(3,2)-box_coords(2,2);
        split_line = round(box_coords(1,2) + (box_size_c_n/2) + (mf*box_size_c_n));

        line([box_coords(1,1),box_coords(2,1)],[box_coords(1,2) box_coords(2,2)],'LineWidth',2,'Color','k'); % Line horizontal low
        line([box_coords(1,1),box_coords(2,1)],[box_coords(3,2) box_coords(4,2)],'LineWidth',2,'Color','k'); % Line horizontal high
        line([box_coords(1,1),box_coords(1,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical left
        line([box_coords(2,1),box_coords(2,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical right
        line([box_coords(1,1),box_coords(2,1)],[split_line split_line],'LineWidth',2,'Color','k'); % Split line

        %set(gca,'xtick',[]);
        %set(gca,'ytick',[]);
        %box off
        xlim([400 1400]);
        ylim([0 1200]);

        title(['Run ' num2str(j)]);

    end
    
    % Save plot
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'Units','centimeters','Position',[10 20 40 15]);
    f = gcf; % f is the handle of the figure you want to export
    if sub(i) < 10
        name = ['Heatmap_trans_Run_Sub0' num2str(sub(i))];
    else
        name = ['Heatmap_trans_Run_Sub' num2str(sub(i))];
    end

    set(gcf,'PaperPositionMode','auto')
    print(f,fullfile(path,name),format,['-r',num2str(rez)])

    close all

end

%% Visualize transformed data per subject (heat maps)
close all

for i = 1:numel(sub)
    
    fprintf('-- Subject %d --\n',sub(i))

    data_s = data_corr(data_corr(:,1)==sub(i),:);
    fix_overall = sum(data_s(:,3));
    
    data_s(data_s(:,4) > monitor_res(1)-area,:) = NaN;
    data_s(data_s(:,5) > monitor_res(2)-area,:) = NaN;

    data_s(data_s(:,4) < area,:) = NaN;
    data_s(data_s(:,5) < area,:) = NaN;

    data_s(any(isnan(data_s),2),:) = [];
    runs = unique(data_s(:,2));
            
    % Heat map template
    heat_map = zeros(monitor_res(2),monitor_res(1));

    % Square area
    for k = 1:size(data_s,1)

        val = data_s(k,6);
        heat_map(data_s(k,5)-area:data_s(k,5)+area,data_s(k,4)-area:data_s(k,4)+area) = heat_map(data_s(k,5)-area:data_s(k,5)+area,data_s(k,4)-area:data_s(k,4)+area) + val;

    end

    % Calculate percent
    heat_map = heat_map/fix_overall*100;

    % Smooth heat map
    g_kernel = fspecial('gaussian',60,60);
    heat_map = filter2(g_kernel,heat_map);

    % Set threshold
    hm_values = sort(reshape(heat_map,[1 size(heat_map,1)*size(heat_map,2)]));
    hm_values(hm_values==0) = [];

    thr_pc = 80; % percent of lower values to be cut off
    vals = unique(hm_values);
    vals_sum = numel(vals);
    val_crit = round(vals_sum*(thr_pc/100));
    cut_off = vals(val_crit);

    % Apply threshold
    heat_map(heat_map < cut_off) = 0;

    %heat_map = heat_map(monitor_res(2)/2-stim_size(2)/2:monitor_res(2)/2+stim_size(2)/2-1,monitor_res(1)/2-stim_size(1)/2:monitor_res(1)/2+stim_size(1)/2-1);

    % Plot
    cm_c = jet;
    cm_c(1,:) = [1 1 1];
    colormap(cm_c);
    imagesc(heat_map);

    hold on
    box_pos_c = [img_center(1)-box_size/2 img_center(2)-box_size/2];
    box_coords = [box_pos_c(1) box_pos_c(2); % top_left
        box_pos_c(1)+box_size box_pos_c(2); % top right
        box_pos_c(1) box_pos_c(2)+box_size; % bottom left
        box_pos_c(1)+box_size box_pos_c(2)+box_size]; % bottom right

    box_coords = round(box_coords * 1.5);
    box_size_c_n = box_coords(3,2)-box_coords(2,2);
    split_line = round(box_coords(1,2) + (box_size_c_n/2) + (mf*box_size_c_n));

    line([box_coords(1,1),box_coords(2,1)],[box_coords(1,2) box_coords(2,2)],'LineWidth',2,'Color','k'); % Line horizontal low
    line([box_coords(1,1),box_coords(2,1)],[box_coords(3,2) box_coords(4,2)],'LineWidth',2,'Color','k'); % Line horizontal high
    line([box_coords(1,1),box_coords(1,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical left
    line([box_coords(2,1),box_coords(2,1)],[box_coords(1,2) box_coords(3,2)],'LineWidth',2,'Color','k'); % Line vertical right
    line([box_coords(1,1),box_coords(2,1)],[split_line split_line],'LineWidth',2,'Color','k'); % Split line

    %set(gca,'xtick',[]);
    %set(gca,'ytick',[]);
    %box off
    xlim([400 1400]);
    ylim([0 1200]);

    % Save plot
    set(gca,'LooseInset',get(gca,'TightInset'));
    set(gcf,'Units','centimeters','Position',[10 20 20 15]);
    f = gcf; % f is the handle of the figure you want to export
    if sub(i) < 10
        name = ['Heatmap_trans_Sub0' num2str(sub(i))];
    else
        name = ['Heatmap_trans_Sub' num2str(sub(i))];
    end

    set(gcf,'PaperPositionMode','auto')
    print(f,fullfile(path,name),format,['-r',num2str(rez)])

    close all

end
