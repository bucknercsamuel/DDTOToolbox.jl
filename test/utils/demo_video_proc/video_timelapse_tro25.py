import cv2
import numpy as np
import ffmpeg
import os
from pdb import set_trace as debug

# NOTE:
# Must install ffmpeg in the following way:
#   > conda install ffmpeg
#   > pip install ffmpeg-python

####################################
### FUNCTION
####################################

# Load video and obtain equally spaced frames from start_time to end_time
def extract_frames(video_path, start_time, end_time, num_frames):
    # Get video duration
    probe = ffmpeg.probe(video_path)
    duration = float(probe['format']['duration'])

    if end_time == -1:
        end_time = duration - 0.1  # Subtract a small value to avoid issue with extracting beyond last frame

    # Calculate time interval between frames
    interval = (end_time - start_time) / (num_frames - 1)

    frames = []
    for i in range(num_frames):
        timestamp = start_time + i * interval
        print(f'Extracted frame {i + 1}/{num_frames} at {timestamp:.2f}s')
        # Extract frame at timestamp
        out, _ = (
            ffmpeg
            .input(video_path, ss=timestamp)
            .output('pipe:', format='image2', vframes=1)
            .run(capture_stdout=True, capture_stderr=True)
        )
        # Convert bytes to numpy array
        frame = np.frombuffer(out, np.uint8)
        
        # Decode image
        ext = os.path.splitext(video_path)[1].lower()
        frame = cv2.imdecode(frame, cv2.IMREAD_COLOR)
        frames.append(frame)
    
    return frames

# Overlay frames such that the moving parts of the frame are fully opaque
def overlay_frames(frames, pixels_diff_threshold=10):
    if not frames:
        raise ValueError("No frames to overlay.")

    # Set first frame as the previous frame
    prev_frame = frames[0].astype(np.float32)

    # Find pixels that have changed between three successive frames
    pixels_diff = []
    for iter in range(1,len(frames) - 2):
        cur_frame = frames[iter]
        next_frame = frames[iter + 1]
        
        # Convert current frame to float32 for precision
        cur_frame = cur_frame.astype(np.float32)

        # Diff: current frame and previous frame
        diff1 = np.abs(cur_frame - prev_frame)
        mask1 = np.any(diff1 > pixels_diff_threshold, axis=-1)
        
        # Diff: next frame and current frame
        diff2 = np.abs(next_frame - cur_frame)
        mask2 = np.any(diff2 > pixels_diff_threshold, axis=-1)
        
        # Combine masks
        mask = np.logical_and(mask1, mask2)

        # Store the mask
        pixels_diff.append(mask)

        # Update previous frame
        prev_frame = cur_frame

    # Create output as overlay of all masks on the first frame
    output_frame = frames[0].copy()
    for frame,mask in zip(frames[1:],pixels_diff):
        # Perform weighted averaging where the mask is true
        alpha = 0.2
        output_frame[mask] = alpha * output_frame[mask] + (1-alpha) * frame[mask]
    
    return output_frame
        
####################################
### MAIN
####################################
if __name__=="__main__": 
    # Parameters for fast drop test
    tag_cam = 'slow_descent'
    t0_cam = 42
    t1_cam = 64
    tag_twin = tag_cam + '_twin'
    t0_twin = 0
    t1_twin = -1
    
    
    
    
    # Parameters, overall
    num_frames = 8  # Number of frames to extract
    num_frames_downsample = 1
    show_timelapse_at_end = False
    separate_bar_width = 5  # Width of the white bar to separate frames
    output_path = 'media/timelapse_' + tag + '.jpg'
    
    # Load raw video
    start_time = 42  # Start time in seconds
    end_time = 64 # End time in seconds (-1 means end of video)
    video_path = 'media/raw/'+tag+'.MOV'  # Path to your video file
    full_video_path = os.path.abspath(video_path)
    frames = extract_frames(full_video_path, start_time, end_time, num_frames)
    timelapse_frame = overlay_frames(frames)
    frames_downsample = frames[::num_frames_downsample]
    
    # Add thin white bar to the right of each downsample frame
    if show_timelapse_at_end:
        offset = 0
    else:
        offset = 1
    for i in range(len(frames_downsample)-offset):
        frames_downsample[i] = np.hstack((frames_downsample[i], np.ones((frames_downsample[i].shape[0], separate_bar_width, 3), dtype=np.uint8) * 255))
    
    # Horizontally concatenate the downsampled frames with the timelapse frame at the end
    if show_timelapse_at_end:
        output_frame_real = np.hstack((np.hstack((frames_downsample)),timelapse_frame))
    else:
        output_frame_real = np.hstack((frames_downsample))
        
    # Load video 2
    start_time = 0  # Start time in seconds
    end_time = -1 # End time in seconds (-1 means end of video)
    video_path = 'media/'+tag+'_twin.mp4'  # Path to your video file
    full_video_path = os.path.abspath(video_path)
    frames = extract_frames(full_video_path, start_time, end_time, num_frames)
    timelapse_frame = overlay_frames(frames)
    frames_downsample = frames[::num_frames_downsample]
    
    # Add thin white bar to the right of each downsample frame
    if show_timelapse_at_end:
        offset = 0
    else:
        offset = 1
    for i in range(len(frames_downsample)-offset):
        frames_downsample[i] = np.hstack((frames_downsample[i], np.ones((frames_downsample[i].shape[0], separate_bar_width, 3), dtype=np.uint8) * 255))
    
    # Add same thin white bar to top of each downsample frame
    for i in range(len(frames_downsample)):
        frames_downsample[i] = np.vstack((np.ones((separate_bar_width, frames_downsample[i].shape[1], 3), dtype=np.uint8) * 255, frames_downsample[i]))
    
    # Horizontally concatenate the downsampled frames with the timelapse frame at the end
    if show_timelapse_at_end:
        output_frame_twin = np.hstack((np.hstack((frames_downsample)),timelapse_frame))
    else:
        output_frame_twin = np.hstack((frames_downsample))
        
    # Concatenate the two video timelapses vertically
    output_frame = np.vstack((output_frame_real, output_frame_twin))
        
    # Save the output frame
    cv2.imwrite(output_path, output_frame)