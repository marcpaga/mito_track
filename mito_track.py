import cv2 as cv
import numpy as np
import tiffcapture as tc
import time
import os
import sys
import copy

# Parameters for Shi-Tomasi corner detection
feature_params = dict(maxCorners = 1000, 
                      qualityLevel = 0.1, 
                      minDistance = 2, 
                      blockSize = 2)
# Parameters for Lucas-Kanade optical flow
lk_params = dict(winSize = (15,15), 
                 maxLevel = 2, 
                 criteria = (cv.TERM_CRITERIA_EPS | cv.TERM_CRITERIA_COUNT, 10, 0.03))
# The video feed is read in as a VideoCapture object
filename = 'dat/flow_wt_4.tif'
name = filename.split('/')[1].split('.')[0]
dots_dir = 'dots/' + name
os.mkdir(dots_dir)
cap = tc.opentiff(filename)
# Variable for color to draw optical flow track
color = (0, 255, 0)
# ret = a boolean return value from getting the frame, first_frame = the first frame in the entire video sequence
ret, first_frame = cap.read()

# Converts frame to grayscale because we only need the luminance channel for detecting edges - less computationally expensive
prev_gray = first_frame

prev_gray = prev_gray.astype(np.float64) / prev_gray.max() # normalize the data to 0 - 1
prev_gray = 255 * prev_gray  # Now scale by 255
prev_gray = prev_gray.astype(np.uint8)


prev = cv.goodFeaturesToTrack(prev_gray, mask = None, **feature_params)
j = 1
new_name = dots_dir + '/new_' + str(j) + '.csv'
np.savetxt(new_name, prev[:,0,:], delimiter = ',')

# Creates an image filled with zero intensities with the same dimensions as the frame - for later drawing purposes
mask = np.zeros_like(first_frame)
cap.release()

cap = tc.opentiff(filename)

while(cap.isOpened()):
    j = j + 1
    
    # ret = a boolean return value from getting the frame, frame = the current frame being projected in the video
    ret, frame = cap.read()
    # Converts each frame to grayscale - we previously only converted the first frame to grayscale

    gray = frame
    gray = gray.astype(np.float64) / gray.max() # normalize the data to 0 - 1
    gray = 255 * gray  # Now scale by 255
    gray = gray.astype(np.uint8)
    # Calculates sparse optical flow by Lucas-Kanade method
    # https://docs.opencv.org/3.0-beta/modules/video/doc/motion_analysis_and_object_tracking.html#calcopticalflowpyrlk
    nextt, status, error = cv.calcOpticalFlowPyrLK(prev_gray, gray, prev, None, **lk_params)
    
    new_name = dots_dir + '/new_' + str(j) + '.csv'
    out_tab = copy.copy(nextt[:,0,:])
    out_tab = np.append(out_tab, status, axis = 1)
    
    np.savetxt(new_name, out_tab, delimiter = ',')
    
    # Selects good feature points for previous position
    good_old = prev[status == 1]
    # Selects good feature points for next position
    good_new = nextt[status == 1]
    
    '''
    if j == 4:
        np.set_printoptions(threshold=sys.maxsize)
        out_tab = copy.copy(nextt[:,0,:])
        print(out_tab.shape)
        print(np.append(out_tab, status, axis = 1))
        exit()
    '''
    
    # Draws the optical flow tracks
    for i, (new, old) in enumerate(zip(good_new, good_old)):
        # Returns a contiguous flattened array as (x, y) coordinates for new point
        a, b = new.ravel()
        # Returns a contiguous flattened array as (x, y) coordinates for old point
        c, d = old.ravel()
        # Draws line between new and old position with green color and 2 thickness
        mask = cv.line(mask, (a, b), (c, d), color, 2)
        # Draws filled circle (thickness of -1) at new position with green color and radius of 3
        frame = cv.circle(frame, (a, b), 3, color, -1)
    # Overlays the optical flow tracks on the original frame
    output = cv.add(frame, mask)
    # Updates previous frame
    prev_gray = gray.copy()
    # Updates previous good feature points
    prev = good_new.reshape(-1, 1, 2)
    
    # Opens a new window and displays the output frame
    # Frames are read by intervals of 10 milliseconds. The programs breaks out of the while loop when the user presses the 'q' key
    cv.imshow('frame',output)
    k = cv.waitKey(30) & 0xff
    time.sleep(0.1)
    if k == 27:
        break
    # Now update the previous frame and previous points
    
# The following frees up resources and closes all windows
cap.release()
cv.destroyAllWindows()