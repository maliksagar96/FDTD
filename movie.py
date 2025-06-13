import cv2
import os

# Frame folder and pattern
frame_folder = 'frames'
frame_pattern = 'frame_{}.png'
start_frame = 0
end_frame = 1000  # inclusive
fps = 30

# Read the first frame to get frame size
first_frame_path = os.path.join(frame_folder, frame_pattern.format(start_frame))
frame = cv2.imread(first_frame_path)
height, width, layers = frame.shape

# Define the codec and create VideoWriter
fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Use 'XVID' for .avi
video = cv2.VideoWriter('movie.mp4', fourcc, fps, (width, height))

for i in range(start_frame, end_frame + 1):
    frame_path = os.path.join(frame_folder, frame_pattern.format(i))
    if os.path.exists(frame_path):
        frame = cv2.imread(frame_path)
        video.write(frame)

video.release()
print("Video saved as movie.mp4")
