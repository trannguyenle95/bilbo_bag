import sys
sys.path.append("/home/david/catkin_ws/src/franka_david/scripts/qdc-manip")

import numpy as np
import cv2
from scipy import ndimage

from franka_real_world.kinect import (KinectClient, compute_coverage, get_cloth_mask,
                                       get_workspace_crop, get_mask_workspace_crop)
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

# depth_width_res = 640
# depth_height_res = 480
# # Availables 1920x1080, 1280x720, 960x540. 6, 15, 30 fps (60 for the latter)
# color_width_res = 1280
# color_height_res = 720
# fps = 30
#
# max_depth = 2.0


gundam_ip = "130.233.123.198"
camera = KinectClient(rotate_img=False, ip=gundam_ip)

print(f"camera intr={camera.color_intr}")
# WS_PC = [128, 569, 122, 750]
# X (128, 750)
# Y (122, 569)

# top_view = cv2.resize(get_workspace_crop(self.top_cam.get_rgbd()[0].copy()),
#                       (256, 256))

while True:
    coverage = compute_coverage(rgb=camera.get_rgbd()[0])
    print(f"Coverage is=[{coverage}]")

    rgb, depth = camera.get_rgbd()
    print(f"shape rgb {rgb.shape}")
    # print(f"RGB shape = {rgb.shape}, depth shape={depth.shape}")

    rgb_to_crop = np.moveaxis(rgb, -1, 0)
    # rotated = ndimage.rotate(rgb, 15)
    cropped_rgb_ws = get_workspace_crop(rgb)
    cropped_rgb = get_mask_workspace_crop(rgb)
    # print(f"shape cropped_rgb {cropped_rgb.shape}")

    cropped_rgb_copy = rgb.copy()
    mask = get_cloth_mask(cropped_rgb)
    print(f"crop dim ={cropped_rgb.shape}")

    # new_exposure = int(input("New exposure="))

    # camera.change_exposure(new_exposure)

    # print(f"Mask shape={mask.shape}")
    # images = np.hstack((rgb[:, :, 0], mask))

    # figure, axes = plt.subplots(1, 2, figsize=(20, 10))
    # axes[0].imshow(rgb)
    # axes[1].imshow(cropped_rgb)
    # plt.show()
    # plt.pause(0.001)
    # plt.clf()

    # Display only the rgb
    # figure, axes = plt.subplots(1, 1, figsize=(10, 10))
    # axes.imshow(rgb)
    # plt.show()
    # plt.pause(0.001)
    # plt.clf()

    # Display both mask and rgb
    figure, axes = plt.subplots(1, 2, figsize=(20, 10))
    axes[0].imshow(mask, cmap='gray', vmin=0, vmax=1.0)
    # print(f"0,0={mask[0,0]}")
    # axes[1].imshow(rgb)
    # figure, axes = plt.subplots(1, 2, figsize=(20, 10))
    # axes[0].imshow(cropped_rgb)
    axes[1].imshow(cropped_rgb)
    # imshow_vir = axes[1].imshow(depth, cmap='viridis', vmin=1.4, vmax=1.5)
    # cbar = figure.colorbar(imshow_vir, ax=axes[1], extend='both')
    plt.show()
    plt.pause(0.001)
    plt.clf()
    plt.close()

    # bottom = (55, 69, 115)
    # top = (237, 233, 227)
    # If we do it with the gray one it is not shown after the mask
    # hsv_image = cv2.inRange(cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV), bottom, top)
    # hsv_image = cv2.inRange(rgb, bottom, top)


    # # Display both mask and rgb
    # figure, axes = plt.subplots(1, 1, figsize=(20, 20))
    # axes.imshow(rgb)
    # # imshow_vir = axes[1].imshow(depth, cmap='viridis', vmin=1.4, vmax=1.5)
    # # cbar = figure.colorbar(imshow_vir, ax=axes[1], extend='both')
    # plt.show()
    # plt.pause(0.001)
    # plt.clf()
    # plt.close()

    # images = np.hstack((cropped_rgb, depth_colormap))
    # inverted_cropped = np.moveaxis(cropped_rgb, -1, 0)
    # print(f"shape inverted crop {inverted_cropped.shape}")
    # Show images

    # scatter_points(info)
    #
    #
    #
    # cv2.namedWindow('RealSense', cv2.WINDOW_AUTOSIZE)
    # depth_colormap = cv2.applyColorMap(cv2.convertScaleAbs(depth, alpha=10.03), cv2.COLORMAP_JET)
    # images = np.hstack((rgb, depth_colormap))
    # imS = cv2.resize(images, (720, 1280))  # Resize image
    # cv2.imshow('Azure Kinect', imS)
    # # cv2.imshow('Azure Kinect', rgb)
    # cv2.waitKey(1)

    # input("Waiting for measuring the coverage again")


# Configure depth and color streams
pipeline = rs.pipeline()
config = rs.config()
# For depth camera available resolutions are: 1024x768=QVGA, 640x480=VGA, 320x240=XGA
# Best are at VGA resolution
depth_width_res = 640
depth_height_res = 480
# Availables 1920x1080, 1280x720, 960x540. 6, 15, 30 fps (60 for the latter)
color_width_res = 1280
color_height_res = 720
fps = 30


# Get device product line for setting a supporting resolution
pipeline_wrapper = rs.pipeline_wrapper(pipeline)
pipeline_profile = config.resolve(pipeline_wrapper)
device = pipeline_profile.get_device()

device_product_line = str(device.get_info(rs.camera_info.product_line))
config.enable_stream(rs.stream.depth, depth_width_res, depth_height_res, rs.format.z16, fps)  # 1024*768 resolution (others are possible)
config.enable_stream(rs.stream.color, depth_width_res, depth_height_res, rs.format.rgb8, fps)



print("Camera configured")
# found_rgb = False
# for s in device.sensors:
#     if s.get_info(rs.camera_info.name) == 'RGB Camera':
#         found_rgb = True
#         break
# if not found_rgb:
#     print("The demo requires Depth camera with Color sensor")
#     exit(0)

print("Enabling stream")
# config.enable_stream(rs.stream.depth, 1024, 768, rs.format.z16, 30)

# Start streaming
pipeline.start(config)
print("Pipeline started")

# Declare depth sensor object and set options
depth_sensor = pipeline_profile.get_device().first_depth_sensor()
# depth_sensor.set_option(rs.option.visual_preset, 5)  # 5 is short range, 3 is low ambient light
# depth_sensor.set_option(rs.option.confidence_threshold, 3)  # 3 is the highest confidence
# depth_sensor.set_option(rs.option.noise_filtering, 6)
# Modified options
depth_sensor.set_option(rs.option.visual_preset, 5)  # 5 is short range, 3 is low ambient light
depth_sensor.set_option(rs.option.confidence_threshold, 1)  # 3 is the highest confidence
# depth_sensor.set_option(rs.option.noise_filtering, 4)
depth_scale = depth_sensor.get_depth_scale()
print(f"Depth Scale is= {depth_scale}")  # 0.0002500000118743628

dist = depth_sensor.set_option(rs.option.min_distance, 1)
dist = depth_sensor.get_option(rs.option.min_distance)
print(f"New min_distance = {dist}")
# depth_sensor.set_option(rs.option.filter_smooth_alpha, 0.17)   # Not available
# depth_sensor.set_option(rs.option.filter_smooth_delta, 48.0)

# Get the depth intrinsics
depth_profile = rs.video_stream_profile(pipeline_profile.get_stream(rs.stream.depth))
depth_intrinsics = depth_profile.get_intrinsics()

# Get the im intrinsics
color_profile = rs.video_stream_profile(pipeline_profile.get_stream(rs.stream.color))
color_intr = color_profile.get_intrinsics()

clipping_distance_in_meters = 1.8  #1 meter
clipping_distance = clipping_distance_in_meters / depth_scale

# Get the depth intrinsics
depth_profile = rs.video_stream_profile(pipeline_profile.get_stream(rs.stream.depth))
depth_intr = depth_profile.get_intrinsics()

# Get the im intrinsics
color_profile = rs.video_stream_profile(pipeline_profile.get_stream(rs.stream.color))
color_intr = color_profile.get_intrinsics()


# Create an align object
# rs.align allows us to perform alignment of depth frames to others frames
# The "align_to" is the stream type to which we plan to align depth frames.
align_to = rs.stream.color
align = rs.align(align_to)

try:
    while True:
        # Get frameset of color and depth
        frames = pipeline.wait_for_frames()

        # Align the depth frame to color frame
        aligned_frames = align.process(frames)

        # Get aligned frames
        aligned_depth_frame = aligned_frames.get_depth_frame()  # aligned_depth_frame is a 640x480 depth image
        color_frame = aligned_frames.get_color_frame()

        # Validate that both frames are valid
        if not aligned_depth_frame or not color_frame:
            continue

        # Convert images to numpy arrays
        depth_image = np.asanyarray(aligned_depth_frame.get_data())
        color_image = np.asanyarray(color_frame.get_data())

        # To get the depth image values in [m]
        depth_scaled = depth_image * depth_scale

        # Remove background - Set pixels further than clipping_distance to grey
        grey_color = 0
        # depth image is 1 channel, color is 3 channels
        depth_image_3d = np.dstack((depth_image, depth_image, depth_image))
        bg_removed = np.where((depth_image_3d > clipping_distance) | (depth_image_3d <= 0), grey_color, color_image)

        # Apply colormap on depth image (image must be converted to 8-bit per pixel first)
        depth_colormap = cv2.applyColorMap(cv2.convertScaleAbs(depth_image, alpha=0.03), cv2.COLORMAP_JET)

        depth_colormap_dim = depth_colormap.shape
        color_colormap_dim = color_image.shape

        # If depth and color resolutions are different, resize color image to match depth image for display
        # if depth_colormap_dim != color_colormap_dim:
        #     resized_color_image = cv2.resize(color_image, dsize=(depth_colormap_dim[1], depth_colormap_dim[0]), interpolation=cv2.INTER_AREA)
        #     images = np.hstack((resized_color_image, depth_colormap))
        # else:
        images = np.hstack((bg_removed, depth_colormap))

        # Show images
        cv2.namedWindow('RealSense', cv2.WINDOW_AUTOSIZE)
        cv2.imshow('RealSense', images)
        cv2.waitKey(1)

finally:
    # Stop streaming
    pipeline.stop()


# # Declare pointcloud object, for calculating pointclouds and texture mappings
# pc = rs.pointcloud()
# # We want the points object to be persistent so we can display the last cloud when a frame drops
# points = rs.points()

# Create pipeline and config stream
# pipeline = rs.pipeline()
# config = rs.config()
# width_res = 1024
# height_res = 768
# fps = 30
# config.enable_stream(rs.stream.depth, width_res, height_res, rs.format.z16, fps) # 1024*768 resolution (others are possible)
# config.enable_stream(rs.stream.color, width_res, height_res, rs.format.rgb8, fps)
#
# profile = config.resolve(pipeline)
# pipeline.start(config)
#
# # Declare depth sensor object and set options
# depth_sensor = profile.get_device().first_depth_sensor()
# depth_sensor.set_option(rs.option.visual_preset, 3)  # 5 is short range, 3 is low ambient light
# depth_sensor.set_option(rs.option.confidence_threshold, 3)  # 3 is the highest confidence
# depth_sensor.set_option(rs.option.noise_filtering, 6)
#
# # Create an align object to match both resolutions
# # The "align_to" is the stream type to which other stream will be aligned.
# align_to = rs.stream.depth
# align = rs.align(align_to)
#
# # Get frames
# frames = pipeline.wait_for_frames() # Wait until a frame is available

# Align the color frame to depth frame
# aligned_frames = align.process(frames)
# depth_frame = aligned_frames.get_depth_frame()
# color_frame = aligned_frames.get_color_frame()