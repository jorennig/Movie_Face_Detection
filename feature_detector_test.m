%% extract images from video
video_list = {'AlexisOhanian'};

for i = 1:numel(video_list)
    
   box_pos = [];
    
    %% EXTRACT FRAMES
    c_dir = pwd;
    name = video_list{i};
    filename = [c_dir '/' name '.mp4'];
    %reading a video file
    mov = VideoReader(filename);
    
    % Defining Output folder as 'snaps'
    opFolder = fullfile(cd, name);
    %if  not existing
    if ~exist(opFolder, 'dir')
        %make directory & execute as indicated in opfolder variable
        mkdir(opFolder);
    end
    
    %getting no of frames
    numFrames = mov.NumberOfFrames;
    
    %setting current status of number of frames written to zero
    numFramesWritten = 0;
    
    %for loop to traverse & process from frame '1' to 'last' frames
    for t = 1:numFrames
        currFrame = read(mov, t);    %reading individual frames
        opBaseFileName = sprintf('%3.3d.png', t);
        opFullFileName = fullfile(opFolder, opBaseFileName);
        imwrite(currFrame, opFullFileName, 'png');   %saving as 'png' file
        %indicating the current progress of the file/frame written
        progIndication = sprintf('Wrote frame %4d of %d.', t, numFrames);
        disp(progIndication);
        numFramesWritten = numFramesWritten + 1;
    end      %end of 'for' loop
    progIndication = sprintf('Wrote %d frames to folder "%s"',numFramesWritten, opFolder);
    disp(progIndication);
    
    main_folder_name = '/Volumes/data/BCM/EyeTracking_Movies/Movie/Movie_Stim/';
    
    %% EXTRACT FEATURES
    for k = 1:numFrames
        
        %figure(k)
        jpgFilename = sprintf('%3.3d.png', k);
        fullFileName = fullfile([main_folder_name '/' name '/'], jpgFilename);
        if exist(fullFileName, 'file')
            I= imread(fullFileName);
        else
            warningMessage = sprintf('Warning: image file does not exist:\n%s', fullFileName);
            uiwait(warndlg(warningMessage));
        end
        
        %% test - find face
        faceDetector = vision.CascadeObjectDetector('MergeThreshold', 20);
        face_bboxes = faceDetector(I); % this gives us the coordinates within the face
        
        %% mouth detector
%         mouthDetector = vision.CascadeObjectDetector('Mouth', 'MergeThreshold', 100);
%         mouth_bboxes = mouthDetector(I); % this gives us the coordinates within the face
        
        %% eyes detector
%         eyesDetector = vision.CascadeObjectDetector('EyePairBig', 'MergeThreshold', 5);
%         eye_bboxes = eyesDetector(I); % this gives us the coordinates within the face
        
        %bound{k} = [face_bboxes; mouth_bboxes; eye_bboxes];
        IFaces{k} = insertObjectAnnotation(I,'rectangle',face_bboxes, 'Face'); hold on
%         IEyes{k} =insertObjectAnnotation(I,'rectangle',eye_bboxes, 'Eyes');
%         IMouth{k} = insertObjectAnnotation(I,'rectangle',mouth_bboxes, 'Mouth');
        
%             subplot(311)
%             imshow(IFaces{k})
%             drawnow
%             title('Face')
%         
        %     subplot(312)
        %     imshow(IMouth{k})
        %     drawnow
        %     title('Mouth')
        %
        %     subplot(313)
        %     imshow(IEyes{k})
        %     drawnow
        %     title('Eyes')
        %
        cFrame_notice = sprintf('Face detected frame %4d of %d.', k, numFrames);
        disp(cFrame_notice);
        
        if isempty(face_bboxes)
            box_pos(k,:) = [NaN, NaN, NaN, NaN];
        else
            box_pos(k,:) = face_bboxes(1,:);
        end
        
    end
    
    save([name '.mat'], 'box_pos');   
end
