%% load in each movie
close all 
clear all

%movies = {'JohnOliver', 'AdamScott','AngelaDuckworth','DonaldGlover','JohnMulaney', 'DanaCarvey', 'TrevorNoah','NickOfferman','OprahWinfrey','RyanReynolds','Retta','Stanley','AmyPoehler', 'TinaFey', 'SalKhan'};
%movies = {'AishaTaylor', 'AubreyPlaza', 'DArcyCarden', 'DaveedDiggs', 'JohnKrasinski', 'MaisieWilliams', 'MargaretBrennan', 'WyattCenac', 'YaraShahidi', 'BobbyFlay'};
%movies = {'GordonRamsey', 'BarackObama', 'KumailNanjiani', 'WayneGretzky', 'JJAbrams', 'MaxGreenfield', 'KristenStewart', 'LillySingh', 'JonStewart', 'AnnDowd', 'MilaKunis', 'EmilyBlunt', 'EllieKemper', 'LucyLiu', 'BenedictCumberbatch', 'Ellen'}
%movies = {'RyanReynolds', 'Stanley', 'WyattCenac', 'DaveedDiggs', 'AishaTaylor'};
%movies = {'AlexisOhanian'};

n_movies = numel(movies);
mf = 0.1; % how far dividing line is from halfway point

for i = 1:n_movies
    
    load(['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Stim/' movies{i} '.mat'])
    n_frames = length(box_pos);
    mean_size = round(nanmean(box_pos(:,3)));
    mean_pos = round(nanmean(box_pos(:,1:2)));
    
    %box_pos = [box_pos(:,1:2) mean_size mean_size];
    
    for j = 1:n_frames         
        %% IDENTIFY SMALL BOX SIZES (probably not correctly on face)
        % -> replace x & y with nearby with gradient (use rounded nums)
        % -> replace width & height with avg box size
        box_size = box_pos(j,3);
        
        if box_size < 350 
            if j== 1
                box_pos(j,:) = [box_pos(j+1,1:2) mean_size mean_size];
            else
                box_pos(j,:) = [box_pos(j-1,1:2) mean_size mean_size];
            end
        end
        
        %% IDENTIFY NANS
        % replace with previous if only 1
        if isnan(box_pos(j,1))
            if ~isnan(box_pos(j+1,1))
                if j == 1
                    box_pos(j,:) = [box_pos(j+1,1:2) mean_size mean_size];
                elseif j < 11 || (j+10 > length(box_pos))
                    box_pos(j,:) = [box_pos(j-1,1:2) mean_size mean_size];
                else
                    nearby = round([nanmean(box_pos(j-10:j+10,1)), nanmean(box_pos(j-10:j+10,2))]);
                    box_pos(j,:) = [nearby mean_size mean_size];
                end
            else
                box_pos(j,:) = [mean_pos mean_size mean_size];
            end
            
        else

        box_pos(j,:) = [box_pos(j,1:2) mean_size mean_size];
        end
       
        %% TEST IF IT'S RIGHT by visualizing
        if j < 10
            img = ['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Stim/' movies{i} '/00' num2str(j) '.png'];
            
        elseif j > 9 && j < 100
            img = ['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Stim/' movies{i} '/0' num2str(j) '.png'];
            
        else
            img = ['/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Stim/'  movies{i} '/' num2str(j) '.png'];
        end

        c_face = box_pos(j,:);
        
        frame_coords = [c_face(1) c_face(2); % top_left
            c_face(1)+c_face(3) c_face(2); % top right
            c_face(1) c_face(2)+c_face(4); % bottom left
            c_face(1)+c_face(3) c_face(2)+c_face(4)]; % bottom right
        
        % Resize frame coords -> monitor size
        % monitor is 1920 x 1080, video is 1280 x 720 --> 1.5 factor
        frame_coords = frame_coords * 1.5;
        
        % Calculate halfway point
        half_dist = (frame_coords(3,2) - frame_coords(1,2))/2;
        halfway_split = frame_coords(1,2) + half_dist + (mf*mean_size);
        
        % Visualize box pos 
        test_box = [frame_coords(1,1) frame_coords(1,2) frame_coords(2,1)-frame_coords(1,1) frame_coords(3,2)-frame_coords(1,2)];
        resized = imresize(imread(img), 1.5);  % img must be resized to show because frame coords now in "monitor" coords
        imshow(resized, [])
        hold on
        rectangle('Position', test_box, 'LineWidth', 2, 'EdgeColor', 'w');
        plot([frame_coords(1,1), frame_coords(2,1)],[halfway_split, halfway_split], 'w', 'LineWidth', 5)
        hold off
        drawnow;
        
    end
   
    save([movies{i} '_adj.mat'],'box_pos')    
end