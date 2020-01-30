# modules
import cv2 as cv
import numpy as np
import tiffcapture as tc
import time
import os
import sys
import copy

## CONFIGURATION ==========
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
filename = 'data/rawdata/flow_wt_1.tif'
output_dir = 'data/processed/'
visualize = True
verbose = True
## CONFIGURATION ==========


name = filename.split('/')
name = name[len(name)-1]
name = name.split('.')[0]
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
    
output_dir_vid = output_dir + name + '/'
if not os.path.isdir(output_dir_vid):
    os.mkdir(output_dir_vid)

   
cap = tc.opentiff(filename)
n_frames = cap.length

## process the first frame
j = 1
if verbose:
    print('Processing frame: ' + str(j) + '/' + str(n_frames), end='\r')    
 
# Variable for color to draw optical flow track
color = (0, 255, 0)
# retrieve first frame
ret, first_frame = cap.read()

# Converts frame to grayscale 
prev_gray = first_frame

# normalize the data to 0 - 1
prev_gray = prev_gray.astype(np.float64) / prev_gray.max() 
# Now scale by 255
prev_gray = 255 * prev_gray  
prev_gray = prev_gray.astype(np.uint8)


prev = cv.goodFeaturesToTrack(prev_gray, mask = None, **feature_params)

# write the point positions
out_file = output_dir_vid + '/points_' + str(j) + '.csv'
np.savetxt(out_file, prev[:,0,:], delimiter = ',')

# Creates an image filled with zero intensities
mask = np.zeros_like(first_frame)
cap.release()

cap = tc.opentiff(filename)

## process all the rest of frames
while(cap.isOpened()):
    j = j + 1
    if verbose:
        print('Processing frame: ' + str(j) + '/' + str(n_frames), end='\r')    

    # get the frame
    ret, frame = cap.read()
    
    # Converts each frame to grayscale 
    gray = frame
    gray = gray.astype(np.float64) / gray.max() # normalize the data to 0 - 1
    gray = 255 * gray  # Now scale by 255
    gray = gray.astype(np.uint8)
    
    # calculate the flow
    nextt, status, error = cv.calcOpticalFlowPyrLK(prev_gray, gray, prev, None, **lk_params)
    
    # write the new point positions
    out_file = output_dir_vid + '/points_' + str(j) + '.csv'
    out_tab = copy.copy(nextt[:,0,:])
    out_tab = np.append(out_tab, status, axis = 1)
    
    np.savetxt(out_file, out_tab, delimiter = ',')
    
    # Selects good feature points for previous position
    good_old = prev[status == 1]
    # Selects good feature points for next position
    good_new = nextt[status == 1]
    
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
    # Frames are read by intervals of 10 milliseconds. 
    if visualize:
        cv.imshow('frame',output)
        k = cv.waitKey(30) & 0xff
        time.sleep(0.1)
    
    if j == n_frames:
        break
    
# The following frees up resources and closes all windows
cap.release()
cv.destroyAllWindows()
