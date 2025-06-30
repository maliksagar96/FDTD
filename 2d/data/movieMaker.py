import cv2
import os
from natsort import natsorted  # to sort files numerically

image_folder = 'frames/Ez/'
video_name = 'ElectricField.mp4'
fps = 30

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
images = natsorted(images)  # sorts like ElectricField_1, ElectricField_2, ...

# Read first image to get dimensions
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

# Define the video writer
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

# Write each frame
for image in images:
    frame = cv2.imread(os.path.join(image_folder, image))
    video.write(frame)

video.release()
print("Video saved:", video_name)
