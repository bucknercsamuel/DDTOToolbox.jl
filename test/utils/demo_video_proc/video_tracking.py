"""
Interactive multi-object tracker (Tier 2 approach).

Workflow:
  1. Opens a video, seeks to `start_time`.
  2. Pops up a window where the user draws N bounding boxes (one per object,
     in the order listed in `object_labels`). SPACE/ENTER confirms each box,
     ESC finishes.
  3. Instantiates one OpenCV CSRT tracker per object, tracks them through
     every frame from `start_time` to `start_time + flight_duration`.
  4. Composes a background image (median-stacked by default for a clean,
     object-free scene) and overlays each object's cubic-spline trajectory
     on top. The spline color is the median BGR color sampled from inside
     each tracker's bounding boxes; a thin white outline behind each spline
     provides contrast against the background.

Requires:
  pip install --upgrade opencv-contrib-python scipy
"""
import os
import cv2
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter


####################################
### FUNCTIONS
####################################

# Return a CSRT tracker instance regardless of OpenCV version layout.
# Modern builds expose cv2.TrackerCSRT_create; legacy contrib exposes cv2.legacy.TrackerCSRT_create.
def create_csrt_tracker():
    if hasattr(cv2, 'legacy') and hasattr(cv2.legacy, 'TrackerCSRT_create'):
        return cv2.legacy.TrackerCSRT_create()
    if hasattr(cv2, 'TrackerCSRT_create'):
        return cv2.TrackerCSRT_create()
    raise RuntimeError(
        "CSRT tracker unavailable. Install opencv-contrib-python: "
        "`pip install --upgrade opencv-contrib-python`"
    )


# Median BGR color across pixels inside all of an object's tracked bboxes.
# An inset reduces background bleed at box edges; bboxes_per_frame entries
# may be None to indicate a frame where the tracker failed.
def median_color_in_bboxes(frames, bboxes_per_frame, inset_frac=0.15):
    samples = []
    for frame, bbox in zip(frames, bboxes_per_frame):
        if bbox is None:
            continue
        x, y, w, h = [int(round(v)) for v in bbox]
        H, W = frame.shape[:2]
        x0, y0 = max(0, x), max(0, y)
        x1, y1 = min(W, x + w), min(H, y + h)
        if x1 <= x0 or y1 <= y0:
            continue
        ix = int(round(inset_frac * (x1 - x0)))
        iy = int(round(inset_frac * (y1 - y0)))
        crop = frame[y0 + iy:max(y0 + iy + 1, y1 - iy), x0 + ix:max(x0 + ix + 1, x1 - ix)]
        if crop.size == 0:
            crop = frame[y0:y1, x0:x1]
        samples.append(crop.reshape(-1, 3))
    if not samples:
        return (255, 255, 255)
    stacked = np.concatenate(samples, axis=0)
    median_bgr = np.median(stacked, axis=0)
    return tuple(int(round(c)) for c in median_bgr)


# Optional Savitzky-Golay smoothing of a 1D sequence, then fit a cubic spline.
# Smoothing tames the per-frame jitter that CSRT inevitably produces.
def smooth_and_spline(ts, vals, smooth_window):
    vals = np.asarray(vals, dtype=np.float64)
    if smooth_window and smooth_window >= 3 and len(vals) > smooth_window:
        w = smooth_window if smooth_window % 2 == 1 else smooth_window + 1
        if w >= 3 and w <= len(vals):
            vals = savgol_filter(vals, w, polyorder=min(3, w - 1))
    return CubicSpline(ts, vals)


# Draw a smooth cubic-spline polyline on `img`. A white outline is drawn first
# (thicker), then the colored line on top (thinner), giving a clean halo
# regardless of background color.
def draw_trajectory(img, ts, xs, ys, color_bgr,
                    n_samples=600, line_thickness=4, outline_thickness=8,
                    smooth_window=7):
    ts = np.asarray(ts, dtype=np.float64)
    if len(ts) < 4:
        pts = np.column_stack([xs, ys]).astype(np.int32)
    else:
        cs_x = smooth_and_spline(ts, xs, smooth_window)
        cs_y = smooth_and_spline(ts, ys, smooth_window)
        t_dense = np.linspace(ts[0], ts[-1], n_samples)
        pts = np.column_stack([cs_x(t_dense), cs_y(t_dense)]).astype(np.int32)

    cv2.polylines(img, [pts], isClosed=False, color=(255, 255, 255),
                  thickness=outline_thickness, lineType=cv2.LINE_AA)
    cv2.polylines(img, [pts], isClosed=False, color=color_bgr,
                  thickness=line_thickness, lineType=cv2.LINE_AA)


# Composite the contents of `frame[bbox]` onto `bg` at the same location,
# using a background-subtraction-derived soft alpha mask so only the moving
# object pixels (not the rectangular bbox surroundings) become visible.
#   bg_full: full-frame median background reference (same shape as `bg`).
#   diff_threshold: per-pixel BGR L2 distance below which we treat the pixel
#       as background (alpha=0).
#   diff_span: alpha ramps from 0 to 1 over this distance above the threshold.
#   alpha_blur_px: gaussian kernel size used to feather the alpha edges.
def composite_object_snapshot(bg, frame, bbox, bg_full,
                              diff_threshold=18.0, diff_span=25.0, alpha_blur_px=7):
    if bbox is None:
        return
    x, y, w, h = [int(round(v)) for v in bbox]
    H, W = bg.shape[:2]
    x0, y0 = max(0, x), max(0, y)
    x1, y1 = min(W, x + w), min(H, y + h)
    if x1 <= x0 or y1 <= y0:
        return

    crop = frame[y0:y1, x0:x1].astype(np.float32)
    bg_crop = bg_full[y0:y1, x0:x1].astype(np.float32)
    diff = np.linalg.norm(crop - bg_crop, axis=-1)
    alpha = np.clip((diff - diff_threshold) / max(diff_span, 1e-6), 0.0, 1.0)

    if alpha_blur_px and alpha_blur_px >= 3:
        k = alpha_blur_px if alpha_blur_px % 2 == 1 else alpha_blur_px + 1
        alpha = cv2.GaussianBlur(alpha, (k, k), 0)

    alpha = alpha[..., None]
    bg[y0:y1, x0:x1] = (alpha * crop + (1.0 - alpha) * bg[y0:y1, x0:x1].astype(np.float32)).astype(np.uint8)


# Median-stack frames to a static "empty arena" background. Sub-samples up to
# `max_samples` evenly across the input to bound memory / time.
def median_background(frames, max_samples=80):
    if len(frames) <= max_samples:
        sel = frames
    else:
        idx = np.linspace(0, len(frames) - 1, max_samples).astype(int)
        sel = [frames[i] for i in idx]
    stack = np.stack(sel, axis=0)
    return np.median(stack, axis=0).astype(np.uint8)


# Enforce strictly increasing timestamps for CubicSpline (in case the video
# decoder returns equal PTS values for consecutive frames).
def enforce_monotonic(ts, eps=1e-6):
    ts = np.asarray(ts, dtype=np.float64).copy()
    for i in range(1, len(ts)):
        if ts[i] <= ts[i - 1]:
            ts[i] = ts[i - 1] + eps
    return ts


####################################
### MAIN
####################################
if __name__ == "__main__":
    # ----- parameters -----
    script_dir = os.path.dirname(os.path.abspath(__file__))
    video_path = os.path.join(script_dir, 'media', 'raw', 'slow_descent.mp4')
    start_time = 3.0          # seconds into the source video
    flight_duration = 33.0    # tracked window length (s)
    background_mode = 'median'  # 'median' | 'first' | 'last'
    output_path = os.path.join(script_dir, 'media', 'tracking_slow_descent.jpg')

    # One ROI is drawn per label, in this order:
    object_labels = ['drone', 'car1', 'car2']

    # Rendering
    line_thickness = 4
    outline_thickness = 9
    spline_samples = 600
    smooth_window = 9         # set to 0 to disable Savitzky-Golay smoothing

    # Stroboscopic object overlay: per object, paste N "ghost" snapshots evenly
    # spaced along its trajectory so the actual object pixels appear at each
    # position. Background subtraction (against the median background) is used
    # to derive a soft alpha mask, so the bbox rectangle never shows.
    show_object_overlays = True
    n_overlay_snapshots = 10  # per object

    # ----- open video & seek to start frame -----
    if not os.path.isfile(video_path):
        raise FileNotFoundError(f"Video file not found: {video_path}")
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        raise RuntimeError(f"Could not open video: {video_path}")

    src_duration = cap.get(cv2.CAP_PROP_FRAME_COUNT) / max(cap.get(cv2.CAP_PROP_FPS), 1e-6)
    end_time = start_time + flight_duration
    if end_time > src_duration:
        print(f"Warning: end_time ({end_time:.2f}s) exceeds video duration "
              f"({src_duration:.2f}s); will track until EOF.")

    cap.set(cv2.CAP_PROP_POS_MSEC, start_time * 1000.0)
    ok, first_frame = cap.read()
    if not ok:
        raise RuntimeError(f"Failed to read frame at t={start_time}s")
    t0_actual = cap.get(cv2.CAP_PROP_POS_MSEC) / 1000.0

    # ----- interactive ROI selection -----
    n_objects = len(object_labels)
    print(f"Draw {n_objects} bounding boxes in this order: {object_labels}")
    print("  Click-drag a box, then press SPACE or ENTER to confirm. "
          "Press ESC after the last box to finish.")
    win = f"Select {n_objects} objects (SPACE=next, ESC=done)"
    bboxes_init = cv2.selectROIs(win, first_frame, showCrosshair=True, fromCenter=False)
    cv2.destroyAllWindows()
    if len(bboxes_init) < n_objects:
        raise RuntimeError(
            f"Expected {n_objects} ROIs, got {len(bboxes_init)}. Aborting."
        )
    bboxes_init = bboxes_init[:n_objects]

    # ----- initialize one tracker per object -----
    trackers = []
    for bbox in bboxes_init:
        tr = create_csrt_tracker()
        tr.init(first_frame, tuple(int(v) for v in bbox))
        trackers.append(tr)

    # ----- track through every frame in the window -----
    timestamps = [t0_actual]
    centers = [[(float(b[0] + b[2] / 2.0), float(b[1] + b[3] / 2.0))] for b in bboxes_init]
    bboxes_per_obj = [[tuple(map(float, b))] for b in bboxes_init]
    frames_used = [first_frame]

    frame_idx = 0
    last_print = 0.0
    while True:
        ok, frame = cap.read()
        if not ok:
            break
        t = cap.get(cv2.CAP_PROP_POS_MSEC) / 1000.0
        if t > end_time:
            break
        frame_idx += 1

        for oi, tr in enumerate(trackers):
            success, bbox = tr.update(frame)
            if success:
                x, y, w, h = bbox
                centers[oi].append((x + w / 2.0, y + h / 2.0))
                bboxes_per_obj[oi].append((x, y, w, h))
            else:
                # Hold last position so the spline stays continuous; mark bbox
                # None so this frame doesn't pollute color sampling.
                centers[oi].append(centers[oi][-1])
                bboxes_per_obj[oi].append(None)

        timestamps.append(t)
        frames_used.append(frame)

        if t - last_print >= 2.0:
            print(f"  Tracked frame {frame_idx} at t={t:.2f}s")
            last_print = t

    cap.release()
    print(f"Total tracked frames: {len(timestamps)}")

    # ----- build background -----
    # We need a clean "empty arena" reference (bg_ref) for the alpha mask used
    # when compositing object snapshots, regardless of which image we actually
    # draw onto (bg).
    needs_median = (background_mode == 'median') or show_object_overlays
    if needs_median:
        print("Computing median background ...")
        bg_ref = median_background(frames_used)
    else:
        bg_ref = first_frame  # not used; only set so name exists

    if background_mode == 'first':
        bg = first_frame.copy()
    elif background_mode == 'last':
        bg = frames_used[-1].copy()
    elif background_mode == 'median':
        bg = bg_ref.copy()
    else:
        raise ValueError(f"Unknown background_mode: {background_mode}")

    # ----- compute per-object color and draw cubic-spline trajectory -----
    ts_arr = enforce_monotonic(timestamps)
    object_colors = []
    for oi, label in enumerate(object_labels):
        xs = [c[0] for c in centers[oi]]
        ys = [c[1] for c in centers[oi]]
        color = median_color_in_bboxes(frames_used, bboxes_per_obj[oi])
        object_colors.append(color)
        print(f"  {label}: median BGR color = {color}")
        draw_trajectory(
            bg, ts_arr, xs, ys, color,
            n_samples=spline_samples,
            line_thickness=line_thickness,
            outline_thickness=outline_thickness,
            smooth_window=smooth_window,
        )

    # ----- composite stroboscopic object snapshots on top of the splines -----
    if show_object_overlays and n_overlay_snapshots > 0 and len(frames_used) >= 2:
        n_snap = min(n_overlay_snapshots, len(frames_used))
        snap_indices = np.linspace(0, len(frames_used) - 1, n_snap).astype(int)
        for oi, label in enumerate(object_labels):
            for fi in snap_indices:
                composite_object_snapshot(
                    bg, frames_used[fi], bboxes_per_obj[oi][fi], bg_ref,
                )
            print(f"  {label}: pasted {n_snap} stroboscopic snapshots")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    cv2.imwrite(output_path, bg)
    print(f"Wrote: {output_path}")
