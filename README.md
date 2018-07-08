# Movie_Face_Detection
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

Fixation_Processing.m\
•	loads in Movie_Fixes.mat & (VideoName).mat
•	for each fix, for each frame: is it on the face or not? is it up or lower?
o	upper = top 65% of box
o	lower = bottom 35%
•	average all frames during fix to get percent time that fix was on the face
•	weight to get where they were looking for entire movie, based on length of each fix
•	save as: MovieFixes_ROI_Tab.mat

ROI_Processing_Movie.m\
asdfa

Transform_Visualize.m\
asdf
