import cv2
import numpy as np
import ffmpeg
import os
from pdb import set_trace as debug

# NOTE:
# Must install ffmpeg in the following way:
#   > conda install ffmpeg
#   > pip install ffmpeg-python

# Parameters
start_time = 59  # Start time in seconds
end_time = 90 # End time in seconds (-1 means end of video)
num_frames = 50  # Number of frames to extract
video_path = 'media/droptest_slow_nonadv.MOV'  # Path to your video file
output_path = 'media/droptest_slow_nonadv_timelapse.jpg'

# Load video and obtain equally spaced frames from start_time to end_time
def extract_frames(video_path, start_time, end_time, num_frames):
    # Get video duration
    probe = ffmpeg.probe(video_path)
    duration = float(probe['format']['duration'])

    if end_time == -1:
        end_time = duration

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
        frame = cv2.imdecode(frame, cv2.IMREAD_COLOR)
        frames.append(frame)
    
    return frames

# Overlay frames such that the moving parts of the frame are fully opaque
pixels_diff = []
pixels_diff_threshold = 30
def overlay_frames(frames):
    if not frames:
        raise ValueError("No frames to overlay.")

    # Set first frame as the previous framee
    prev_frame = frames[0].astype(np.float32)

    # Find pixels that have changed between three successive frames
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
        

full_video_path = os.path.abspath(video_path)
frames = extract_frames(full_video_path, start_time, end_time, num_frames)
output_frame = overlay_frames(frames)

# Save the output frame
cv2.imwrite(output_path, output_frame)