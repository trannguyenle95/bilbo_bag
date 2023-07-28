
import requests
import pickle
import numpy as np
import cv2
from scipy import ndimage
from threading import Thread
import skimage.morphology as morph

import time
# import socket
# import threading
# See https://github.com/columbia-ai-robotics/PyKinect


from deps.flingbot.environment.utils import get_largest_component
from franka import WS_PC, WS_PC_MASK


def get_largest_component(arr):
    # label connected components for mask
    #TODO Fix, if there are no components we will get None
    labeled_arr, num_components = morph.label(arr, return_num=True, background=0)
    masks = [(i, (labeled_arr == i).astype(np.uint8)) for i in range(0, num_components)]
    masks.append((len(masks), 1-(np.sum(mask for i, mask in masks) != 0)))
    sorted_volumes = sorted(masks, key=lambda item: np.count_nonzero(item[1]), reverse=True)
    for i, mask in sorted_volumes:
        if arr[mask == 1].sum() == 0:
            continue
        return mask


def setup_thread(target):
    thread = Thread(target=target)
    thread.daemon = True
    thread.start()
    return thread


def get_workspace_crop(img):
    # print(f"Image shape = {img.shape}")
    retval = img[WS_PC[0]:WS_PC[1], WS_PC[2]:WS_PC[3], ...]
    return retval


def get_mask_workspace_crop(img):
    # print(f"Image shape = {img.shape}")
    retval = img[WS_PC_MASK[0]:WS_PC_MASK[1], WS_PC_MASK[2]:WS_PC_MASK[3], ...]
    return retval


def pix_to_3d_position(
        x, y, depth_image, cam_intr, cam_extr, cam_depth_scale):
    # Get click point in camera coordinates
    click_z = depth_image[y, x]* cam_depth_scale
    click_x = (x-cam_intr[0, 2]) * \
        click_z/cam_intr[0, 0]
    click_y = (y-cam_intr[1, 2]) * \
        click_z/cam_intr[1, 1]
    if click_z == 0:
        raise InvalidDepthException
    click_point = np.asarray([click_x, click_y, click_z])
    click_point = np.append(click_point, 1.0).reshape(4, 1)

    # Convert camera coordinates to robot coordinates
    target_position = np.dot(cam_extr, click_point)
    target_position = target_position[0:3, 0]
    return target_position


def get_cloth_mask(rgb):
    h, w, c = rgb.shape
    # Rotated image 1027x1423
    # USE WS_PC that does not use the manipulation one
    if h == 720 and w == 1280 or h == 1027 and w == 1423 or h == 1080 and w == 1920:
        # print("Shape of image")
        rgb[:WS_PC_MASK[0], ...] = 1
        rgb[WS_PC_MASK[1]:, ...] = 1
        rgb[:, :WS_PC_MASK[2], :] = 1
        rgb[:, WS_PC_MASK[3]:, :] = 1
    """
    Segments out black backgrounds
    """

    # Checkered rag  ---- old foam
    # bottom = (68, 23, 89)
    # top = (255, 255, 255)

    # HSV 0-1, 0-1, 0.592-1.0 , 0.592=150.96

    # RGB Gray foam
    # bottom = (0, 0, 170)  # Should be 151, 170 removes a bit
    # bottom = (0, 0, 150)  # Should be 151, 170 removes a bit
    # top = (255, 255, 255)
    # bottom = (110, 0, 78)  # 155, 0 ,78 // 170 with all light
    # huge reflection. 110 keeps some gray

    # top = (255, 255, 255)
    # mask = cv2.inRange(rgb, bottom, top)
    # HSV

    # 20th February 2020 - new lab lighting
    # bottom = (143, 137, 0)  # Should be 151, 170 removes a bit
    # top = (255, 255, 255)
    # mask = cv2.inRange(rgb, bottom, top)

    bottom = (0.0, 0.0, 0.562*255)  # Should be 151, 170 removes a bit
    top = (255, 255, 255)
    mask = cv2.inRange(cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV), bottom, top)

    mask = (mask == 0).astype(np.uint8)
    mask = (1 - mask)  # Need to invert the mask
    if mask.shape[0] != mask.shape[1]:
        mask[:, :int(mask.shape[1]*0.2)] = 0
        mask[:, -int(mask.shape[1]*0.2):] = 0
    return get_largest_component(mask).astype(np.uint8)


def compute_coverage(rgb):
    mask = get_cloth_mask(rgb=rgb)
    return np.count_nonzero(mask) / (mask.shape[0] * mask.shape[1])


class InvalidDepthException(Exception):
    def __init__(self):
        super().__init__('Invalid Depth Point')


class KinectClient:
    def __init__(self, rotate_img, ip='localhost', port=8080):
        self.ip = ip
        self.port = port
        self.rotation_angle = 15
        self.rotate_img = rotate_img
        # self.rotated_h = 1080
        # self.rotated_w = 1920
        self.rotated_h = 720
        self.rotated_w = 1280

    @property
    def color_intr(self):
        color_intr = self.get_intr()
        color_intr[0, 2] = self.rotated_w/2
        color_intr[1, 2] = self.rotated_h/2
        return color_intr
        # if self.rotate_img:
        #     color_intr = self.get_intr()
        #     color_intr[0, 2] = self.rotated_w/2
        #     color_intr[1, 2] = self.rotated_h/2
        #     return color_intr
        # else:
        # return self.get_intr()

    def get_intr(self):
        return pickle.loads(requests.get(f'http://{self.ip}:{self.port}/intr').content)

    def change_exposure(self, exposure):

        if exposure < 0:
            print("Setting auto-exposure")
            requests.get(f'http://{self.ip}:{self.port}/exposure?auto=True&shutter_us={8330}')

        else:
            end_exposure = min(10300, exposure)
            print(f"Setting exposure to ={end_exposure}")

            requests.get(f'http://{self.ip}:{self.port}/exposure?auto=False&shutter_us={end_exposure}')

    def change_brightness(self, brightness):

        end_brightness = min(255, brightness)
        end_brightness = max(0, end_brightness)
        print(f"Setting brightness to ={brightness}")

        requests.get(f'http://{self.ip}:{self.port}/exposure?brightness={end_brightness}')

    def get_rgbd(self, repeats=2):
        time.sleep(0.05)
        data = pickle.loads(requests.get(
            f'http://{self.ip}:{self.port}/pickle/{repeats}').content)
        if self.rotate_img:
            color_rotated = ndimage.rotate(data['color_img'], self.rotation_angle)
            depth_rotated = ndimage.rotate(data['depth_img'], self.rotation_angle)
            return color_rotated, depth_rotated
        else:

            # define the alpha and beta
            # alpha = 1.2  # Contrast control
            # beta = 10.0  # Brightness control

            # call convertScaleAbs function
            color_img = data['color_img']
            # adjusted = cv2.convertScaleAbs(color_img, alpha=alpha, beta=beta)
            # adjusted = cv2.addWeighted( color_img, alpha, color_img, 0, beta)
            #
            return color_img, data['depth_img']
            # return adjusted, data['depth_img']


