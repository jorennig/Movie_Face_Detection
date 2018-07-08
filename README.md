### Movie_Face_Detection
Detects faces in videos, creates faces ROIs per frame, corrects outliers, summarizes data over movie in variables and heat maps

feature_detector_test.m\
•	Extract frames from mp4s & saves as png
•	Get face position in each frame using computer vision toolbox Cascade Object Detector
o	box_pos = [top left x, top left y, width, height]
o	No face detected, box_pos = NaNs
•	Save as (PersonName).mat

face_box_adjustment.m\
•	loads in matrices with box pos  from(PersonName).mat
•	id frames with small face box sizes size or NaN and make pos = nearby box positions
•	change size for all boxes to the mean size for video
•	visualize where the box is on the face to check
o	coords need to be multiplied by 1.5 to make them compatible to the monitor (which is where the ET coordinates are going to come from)
o	The image is also sized up here to simulate the monitor
•	Save as (PersonName)_adj.mat

EyeMovement_Processing.m\
•	ascii ET data --> list of fixations with x, y, duration, start, end time
•	save as Movie_Fixes.mat

ROI_Processing_Movie.m\
•	loads in Movie_Fixes.mat & (VideoName).mat
•	for each fix, for each frame: is it on the face or not? is it up or lower?
o	upper = top 60% of box
o	lower = bottom 40%
•	average all frames during fix to get percent time that fix was on the face
•	weight to get where they were looking for entire movie, based on length of each fix
•	save as: MovieFixes_ROI_Tab.mat

Transform_Visualize.m\
Shows all fixations per movie stimulus for each subject.\
Finds all frames that were shown during a fixation time, adjusts their position relative to average frame, centers the data and visualizes individual transfomed fixations per movie.\
Summarizes transformed/normalized viewing behavior across movies per run/subject.\
