import cv2
import numpy as np
import ffmpeg
import os
from PIL import Image, ImageDraw, ImageFont
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
    if not os.path.isfile(video_path):
        raise FileNotFoundError(f"Video file not found: {video_path}")

    # Get video duration. Re-raise with ffprobe's stderr so the failure is diagnosable.
    try:
        probe = ffmpeg.probe(video_path)
    except ffmpeg.Error as e:
        stderr = e.stderr.decode('utf-8', errors='replace') if e.stderr else '(no stderr)'
        raise RuntimeError(f"ffprobe failed for {video_path}:\n{stderr}") from e
    duration = float(probe['format']['duration'])

    if end_time == -1:
        end_time = duration - 0.1  # Subtract a small value to avoid issue with extracting beyond last frame

    # Calculate time interval between frames
    interval = (end_time - start_time) / (num_frames - 1)

    frames = []
    for i in range(num_frames):
        timestamp = start_time + i * interval
        print(f'Extracted frame {i + 1}/{num_frames} at {timestamp:.2f}s')
        # Extract frame at timestamp. Use image2pipe + mjpeg because the bare
        # image2 muxer is a file-writer and can silently emit zero bytes to a pipe.
        try:
            out, err = (
                ffmpeg
                .input(video_path, ss=timestamp)
                .output('pipe:', format='image2pipe', vcodec='mjpeg', vframes=1)
                .run(capture_stdout=True, capture_stderr=True)
            )
        except ffmpeg.Error as e:
            stderr = e.stderr.decode('utf-8', errors='replace') if e.stderr else '(no stderr)'
            raise RuntimeError(
                f"ffmpeg failed extracting frame at {timestamp:.2f}s from {video_path}:\n{stderr}"
            ) from e

        if not out:
            stderr = err.decode('utf-8', errors='replace') if err else '(no stderr)'
            raise RuntimeError(
                f"ffmpeg returned empty buffer at {timestamp:.2f}s from {video_path} "
                f"(video duration={duration:.2f}s). ffmpeg stderr:\n{stderr}"
            )

        frame_arr = np.frombuffer(out, np.uint8)
        frame = cv2.imdecode(frame_arr, cv2.IMREAD_COLOR)
        if frame is None:
            raise RuntimeError(f"cv2.imdecode failed at {timestamp:.2f}s from {video_path}")
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

# Resolve serif italic + roman font paths once. On Windows these ship with the
# OS; on other platforms fall back to PIL's default bitmap font.
def _resolve_serif_fonts():
    candidates_italic = [
        r"C:\Windows\Fonts\timesi.ttf",                 # Times New Roman Italic (Windows)
        "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
        "/Library/Fonts/Times New Roman Italic.ttf",
    ]
    candidates_roman = [
        r"C:\Windows\Fonts\times.ttf",                  # Times New Roman (Windows)
        "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
        "/Library/Fonts/Times New Roman.ttf",
    ]
    italic = next((p for p in candidates_italic if os.path.isfile(p)), None)
    roman = next((p for p in candidates_roman if os.path.isfile(p)), None)
    return italic, roman

_ITALIC_FONT_PATH, _ROMAN_FONT_PATH = _resolve_serif_fonts()

# Render a LaTeX-style time label (italic math + roman units) to a BGR uint8
# image of size (height_px, width_px) with a white background. Mimics
# "$t = X.X$ s": italic serif for "t = X.X", roman serif for the trailing " s".
# border_px: width (in pixels) of an optional black border drawn just inside
# the label box. Set to 0 to disable.
def render_latex_label(math_part, unit_part, width_px, height_px, border_px=0):
    img = Image.new('RGB', (width_px, height_px), color=(255, 255, 255))
    draw = ImageDraw.Draw(img)

    # Auto-fit font size so the combined string fits within ~92% of the box.
    target_h = int(height_px * 0.72)
    if _ITALIC_FONT_PATH is not None and _ROMAN_FONT_PATH is not None:
        font_italic = ImageFont.truetype(_ITALIC_FONT_PATH, target_h)
        font_roman = ImageFont.truetype(_ROMAN_FONT_PATH, target_h)
    else:
        # Fallback (rare; will look less LaTeX-y but still legible)
        font_italic = ImageFont.load_default()
        font_roman = font_italic

    def measure(font, text):
        bbox = font.getbbox(text)
        return bbox[2] - bbox[0], bbox[3] - bbox[1]

    w_math, h_math = measure(font_italic, math_part)
    w_unit, h_unit = measure(font_roman, unit_part)
    total_w = w_math + w_unit
    max_w = int(width_px * 0.92)

    # Shrink if too wide
    if total_w > max_w:
        scale = max_w / total_w
        target_h = max(8, int(target_h * scale))
        font_italic = ImageFont.truetype(_ITALIC_FONT_PATH, target_h) if _ITALIC_FONT_PATH else font_italic
        font_roman = ImageFont.truetype(_ROMAN_FONT_PATH, target_h) if _ROMAN_FONT_PATH else font_roman
        w_math, h_math = measure(font_italic, math_part)
        w_unit, h_unit = measure(font_roman, unit_part)
        total_w = w_math + w_unit

    # Use ascent for stable vertical alignment between italic and roman glyphs.
    ascent_italic, descent_italic = font_italic.getmetrics() if hasattr(font_italic, 'getmetrics') else (h_math, 0)
    ascent_roman, descent_roman = font_roman.getmetrics() if hasattr(font_roman, 'getmetrics') else (h_unit, 0)
    line_h = max(ascent_italic + descent_italic, ascent_roman + descent_roman)

    x0 = (width_px - total_w) // 2
    y0 = (height_px - line_h) // 2
    # Account for each font's internal bbox origin (getbbox()[0] is left bearing)
    lb_math = font_italic.getbbox(math_part)[0]
    lb_unit = font_roman.getbbox(unit_part)[0]
    draw.text((x0 - lb_math, y0), math_part, font=font_italic, fill=(0, 0, 0))
    draw.text((x0 + w_math - lb_unit, y0), unit_part, font=font_roman, fill=(0, 0, 0))

    # Optional black border drawn inside the label box. Pillow's rectangle width
    # grows symmetrically about the edge, so inset by half so it stays inside.
    if border_px > 0:
        b = max(1, int(border_px))
        inset = b // 2
        draw.rectangle(
            [inset, inset, width_px - 1 - inset, height_px - 1 - inset],
            outline=(0, 0, 0), width=b,
        )

    rgb = np.array(img)
    bgr = cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR)
    return bgr

# Overlay a LaTeX-formatted time indicator on each frame and return new list.
# times_s: iterable of times (s) matching len(frames). height_frac/width_frac
# size the label box relative to each frame; y_frac is the vertical center of
# the label as a fraction of frame height (0.85 ≈ near the bottom).
def add_time_labels(frames, times_s, height_frac=0.1, width_frac=0.98, y_frac=2.0, border_px=0):
    out = []
    for frame, t in zip(frames, times_s):
        H, W = frame.shape[:2]
        h = int(round(H * height_frac))
        w = int(round(W * width_frac))
        x0 = (W - w) // 2
        y0 = int(round(H * y_frac)) - h // 2
        # Clamp to frame bounds
        x0 = max(0, min(x0, W - w))
        y0 = max(0, min(y0, H - h))
        math_part = f"t = {t:.1f}"
        unit_part = " s"
        label_img = render_latex_label(math_part, unit_part, w, h, border_px=border_px)
        f = frame.copy()
        f[y0:y0 + h, x0:x0 + w] = label_img
        out.append(f)
    return out

####################################
### MAIN
####################################
if __name__=="__main__":
    # Parameters for fast drop test
    # Set tag_cam to None to skip the raw camera video (e.g. simulation-only runs)
    tag_cam = 'slow_descent'
    t0_cam = 3
    tag_twin = 'slow_descent_twin'
    t0_twin = 3
    flight_duration = 33  # Real-time flight duration shared by both videos (seconds)

    # Parameters, overall
    num_frames = 8  # Number of frames to extract
    num_frames_downsample = 1
    show_timelapse_at_end = False
    show_time_labels = True  # Overlay LaTeX '$t = X.X$ s' on the bottom row of frames
    label_border_px = 3      # Width (px) of optional black border around each label; set to 0 to disable
    separate_bar_width = 5  # Width of the white bar to separate frames

    # Relative flight times for each displayed downsample frame (used for time labels).
    # Labels run from t=0 at the first frame to t=flight_duration at the last,
    # mirroring how t0_cam / t0_twin are aligned to the start of the flight.
    times_displayed = np.linspace(0, flight_duration, num_frames)[::num_frames_downsample]
    output_tag = tag_cam if tag_cam is not None else tag_twin
    # Anchor all relative paths to this script's directory so it works regardless of CWD.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(script_dir, 'media', 'timelapse_' + output_tag + '.jpg')

    output_frame_real = None

    # Load raw video (optional)
    if tag_cam is not None:
        t1_cam = t0_cam + flight_duration
        full_video_path = os.path.join(script_dir, 'media', 'raw', tag_cam + '.mp4')
        frames = extract_frames(full_video_path, t0_cam, t1_cam, num_frames)
        timelapse_frame = overlay_frames(frames)
        frames_downsample = frames[::num_frames_downsample]

        # Add thin white bar to the right of each downsample frame
        if show_timelapse_at_end:
            offset = 0
        else:
            offset = 1
        for i in range(len(frames_downsample) - offset):
            frames_downsample[i] = np.hstack((frames_downsample[i], np.ones((frames_downsample[i].shape[0], separate_bar_width, 3), dtype=np.uint8) * 255))

        # Horizontally concatenate the downsampled frames with the timelapse frame at the end
        if show_timelapse_at_end:
            output_frame_real = np.hstack((np.hstack((frames_downsample)), timelapse_frame))
        else:
            output_frame_real = np.hstack((frames_downsample))

    # Load twin video
    t1_twin = t0_twin + flight_duration
    full_video_path = os.path.join(script_dir, 'media', 'raw', tag_twin + '.mp4')
    frames = extract_frames(full_video_path, t0_twin, t1_twin, num_frames)
    timelapse_frame = overlay_frames(frames)
    frames_downsample = frames[::num_frames_downsample]

    # Overlay LaTeX time labels on the bottom (twin) row when enabled — the twin
    # row is always the bottom row in the final composite, regardless of cam presence.
    if show_time_labels:
        frames_downsample = add_time_labels(frames_downsample, times_displayed, border_px=label_border_px)

    # Add thin white bar to the right of each downsample frame
    if show_timelapse_at_end:
        offset = 0
    else:
        offset = 1
    for i in range(len(frames_downsample) - offset):
        frames_downsample[i] = np.hstack((frames_downsample[i], np.ones((frames_downsample[i].shape[0], separate_bar_width, 3), dtype=np.uint8) * 255))

    # Add same thin white bar to top of each downsample frame and the timelapse frame
    # (only needed when stacked under the cam row). Applying it to the timelapse too keeps
    # heights consistent for the show_timelapse_at_end hstack below.
    if output_frame_real is not None:
        for i in range(len(frames_downsample)):
            frames_downsample[i] = np.vstack((np.ones((separate_bar_width, frames_downsample[i].shape[1], 3), dtype=np.uint8) * 255, frames_downsample[i]))
        timelapse_frame = np.vstack((np.ones((separate_bar_width, timelapse_frame.shape[1], 3), dtype=np.uint8) * 255, timelapse_frame))

    # Horizontally concatenate the downsampled frames with the timelapse frame at the end
    if show_timelapse_at_end:
        output_frame_twin = np.hstack((np.hstack((frames_downsample)), timelapse_frame))
    else:
        output_frame_twin = np.hstack((frames_downsample))

    # Concatenate the two video timelapses vertically, or output only the twin if no cam video
    if output_frame_real is not None:
        output_frame = np.vstack((output_frame_real, output_frame_twin))
    else:
        output_frame = output_frame_twin

    # Save the output frame
    cv2.imwrite(output_path, output_frame)