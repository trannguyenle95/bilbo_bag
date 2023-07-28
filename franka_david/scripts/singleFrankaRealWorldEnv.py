import os
import cv2
import h5py
import numpy as np
import torch

from filelock import FileLock
from copy import deepcopy
from time import time, sleep, strftime
from PIL import Image
import matplotlib
from matplotlib import pyplot as plt
import hashlib
import imageio
import tqdm

matplotlib.use('TkAgg')


from deps.flingbot.learning.nets import prepare_image
from deps.flingbot.environment.utils import (preprocess_obs, add_text_to_image,
                                             pixels_to_3d_positions, visualize_action)
from deps.flingbot.environment.tasks import Task


from envs.SequenceMemory import SequentialMemory
# from envs.singleArmEnv import SingleArmEnv  # Should maybe use the episodicSingleArmEnv
from agents.agent_utils import (pixels_to_3d_positions_filter_place, get_closest_2d,
                                pixel_to_3d, get_pretransformed_pixels)
from envs.episodEvalSingleArmEnv import evalSingleArmEnv
from envs.plot_utils import visualize_action_w_place


from franka_real_world.kinect import (KinectClient, get_workspace_crop, pix_to_3d_position,
                                      compute_coverage, get_cloth_mask, setup_thread,
                                      get_mask_workspace_crop)

from franka_real_world.franka import (Franka, DEFAULT_ORN, WS_PC, bound_grasp_pos)
# from real_world.realur5_utils import setup_thread
from real_world.real_setup import CLOTHS_DATASET, ask_for_new_cloth
from real_world.realsense_l515 import RealSense, get_video_crop

from envs.primitives import (get_action_params, check_action_reachability,
                             check_arm_reachability, get_xyz_coordinates_given_mid_theta_ref)


class GraspFailException(Exception):
    def __init__(self):
        super().__init__('Grasp failed due to real world')


class FrankaRealWorldEnv(evalSingleArmEnv):
    def __init__(self, randomise_rotation, replace_background=True, **kwargs):
        self.replace_background = replace_background

        # TODO Modify this function to ask for cloth grasped and move to the top, then to the bottom and drop
        def randomize_cloth():
            cloth_type = ask_for_new_cloth()
            self.franka.out_of_the_way()
            rgb_image = self.top_cam.get_rgbd()[0]
            mask_initial = get_cloth_mask(rgb_image.copy())
            figure, axes = plt.subplots(1, 2, figsize=(20, 10))
            np.save('rgb_image_small_towel.npy', rgb_image)
            axes[0].imshow(rgb_image)
            axes[1].imshow(mask_initial)
            # axes[1].matshow(mask)
            plt.show()
            plt.pause(0.001)
            print("Creating new task")
            return Task(
                name=f'{cloth_type}' + strftime("%Y-%m-%d_%H-%M-%S"),
                flatten_area=CLOTHS_DATASET[cloth_type]['flatten_area'],
                initial_coverage=compute_coverage(rgb=rgb_image),
                task_difficulty='hard',
                cloth_mass=CLOTHS_DATASET[cloth_type]['mass'],
                cloth_size=CLOTHS_DATASET[cloth_type]['cloth_size'],
            )

        super().__init__(get_task_fn=randomize_cloth, parallelize_prepare_image=False,
                         **kwargs)

        # state variables to handle recorder
        self.recording = False
        self.recording_daemon = None
        self.recording_daemon_side = None
        # np.random.seed(int(time()))
        self.number_of_restarts = 0
        if randomise_rotation:
            self.random_rotation = np.array([0.0])
            self.restart_rotations = np.array([(i * 15) for i in range(13)])  # From 0 to 180
        else:
            self.restart_rotations = np.array([(i * 15) for i in range(13)])  # From 0 to 180

        # Modify the intrinsics matrix
        self.intrinsics_matrix = self.top_cam.color_intr

        self.depth_limit = 1.425

        # self.ws_pc = WS_PC
        self.action_handlers = {
            'dynamic': self.dynamic_primitive,
            'drag': self.pick_and_drag_primitive,
            'place': self.pick_and_place_primitive,
            'sequential-dynamic': self.dynamic_primitive,
            'sequential-drag': self.pick_and_drag_primitive,
            'sequential-place': self.pick_and_place_primitive,
            'sequential-place-lift': self.pick_and_place_primitive
        }

    def setup_env(self):
        # TODO We will use only one camera
        # used for getting obs
        gundam_ip = "130.233.123.198"

        self.top_cam = KinectClient(rotate_img=False, ip=gundam_ip)

        # used for recording visualizations
        # if you have a third camera/webcam,
        # you can setup a different camera here
        # with a better view of both arms
        # Availables 1920x1080, 1280x720, 960x540. 6, 15, 30 fps (60 for the latter)
        depth_width_res = 640
        depth_height_res = 480
        # Availables 1920x1080, 1280x720, 960x540. 6, 15, 30 fps (60 for the latter)
        color_width_res = 1280
        color_height_res = 720
        #  Color 960x540 @ 30Hz BGR8
        fps = 30

        max_depth = 2.0

        self.setup_cam = RealSense(max_depth, im_h=color_height_res, im_w=color_width_res,
                                   depth_h=depth_height_res, depth_w=depth_width_res, fps=fps)

        # self.setup_cam = None
        self.franka = Franka(init_node=True)
        # self.franka.send_home()
        # self.franka.out_of_the_way()

        self.root_dir = "/home/david/catkin_ws/src/franka_david/scripts/qdc-manip/franka_real_world/"
        cam_pose_txt_file = "top_down_left_ur5_cam_pose.txt"
        # cam_pose_txt_file = "ros_cam_pose.txt"
        cam_depth_file = "camera_depth_scale.txt"
        self.top_cam_left_ur5_pose = np.loadtxt(self.root_dir+cam_pose_txt_file)
        self.pose_matrix = np.loadtxt(self.root_dir+cam_pose_txt_file)
        # self.cam_depth_scale = self.top_cam.depth_scale
        self.cam_depth_scale = np.loadtxt(self.root_dir+cam_depth_file)

    def get_cloth_mask(self, rgb=None):
        if rgb is None:
            rgb = self.top_cam.get_rgbd()[0]
        return get_cloth_mask(rgb)

    def preaction(self):
        self.preaction_mask = self.get_cloth_mask()

    def compute_iou(self):
        mask = self.get_cloth_mask()
        intersection = np.logical_and(mask, self.preaction_mask).sum()
        union = np.logical_or(mask, self.preaction_mask).sum()
        return intersection/union

    def postaction(self):
        curr_coverage = self.compute_coverage()
        self.check_terminal_state(curr_coverage)

    def ask_termination(self):
        if self.current_timestep >= 10:
            print(f"Checking termination. Current episode is={self.current_timestep}")
            terminate = input("Have 10 episodes terminated without grasp errors?, use 0/1")
            if terminate == '0':
                return False
            elif terminate == '1':
                return True
        else:
            return False

    def step(self, value_maps):
        # NOTE: negative current_timestep's
        # act as error codes
        print(f'Step {self.current_timestep}')
        try:
            # if 'sequential' in self.action_primitive[0]:
            #     list_pick_value_maps = []
            #     for action in value_maps:
            #         _, _, pick_value_maps, _, _ = action[0]
            #         max_q_value = pick_value_maps.cpu()[0, :, :].numpy().max()
            #         list_pick_value_maps.append(max_q_value)
            #     sorted_q_idx = np.flip(np.array(list_pick_value_maps).argsort())
            #     for i in range(sorted_q_idx.shape[0]):
            #         selected_next = sorted_q_idx[i]
            #         rotation_angle = self.restart_rotations[i]
            #         # Modify the rotations list so that we use the correct rotation for the check pick place action
            #         self.rotations = [rotation_angle]
            #         next_value_map = value_maps[selected_next][0]
            #         retval = super().step(next_value_map)  # Performs the singleArmEnv
            #         if self.number_of_restarts == 0:
            #             break
            # else:
            retval = super().step(value_maps)  # Performs the singleArmEnv
            self.episode_memory.add_value(key='raw_observation', value=self.raw_pretransform_rgb)
            self.episode_memory.add_value(key='failed_grasp', value=0)
            self.episode_memory.add_value(key='timed_out', value=0)
            self.episode_memory.add_value(key='cloth_stuck', value=0)
            self.episode_memory.add_value(key='task_name', value=self.current_task.name)
            return retval

        except GraspFailException as e:
            # action failed in real world
            print('\t[ERROR]', e)
            self.current_timestep -= 1

            if self.dump_visualizations:
                self.stop_recording()

            self.episode_memory.data['failed_grasp'] = [1] * len(self.episode_memory)

            obs, num_components = self.get_obs()  # Add next observation transformed

            transformed_next_obs = prepare_image(obs, self.get_transformations(), self.obs_dim,
                                                 parallelize=self.parallelize_prepare_image)

            if self.grayscale_only:  # Let's transform to grayscale the next observation
                gray_next_obs = self.rgbd_to_grayscale(transformed_next_obs[0])
                self.episode_memory.add_value(key='next_observations', value=gray_next_obs[0])
            else:
                self.episode_memory.add_value(key='next_observations', value=transformed_next_obs[0])

            if self.grayscale_only:
                return gray_next_obs, self.ray_handle
            else:
                return self.transformed_obs, self.ray_handle

    def ask_action_valid(self):
        valid_action = input("Is the input valid?, use 0/1")
        if valid_action == '0':
            return False
        elif valid_action == '1':
            return True
        else:
            valid_action = input("Is the input valid?, use 0/1")
            if valid_action == '0':
                return False
            elif valid_action == '1':
                return True
            else:
                raise ValueError("Wrong input by the user, should have written either 0 or 1")

    def visualize_and_ask_if_valid(self, retval):

        action_vis = retval['get_action_visualization_fn']()
        # action_vis_array = np.array(action_vis)
        # img_action_visualisation = np.swapaxes(np.swapaxes(action_vis_array, -1, 0), 1, 2)
        im = Image.fromarray(action_vis)
        im.show()

        # figure, axes = plt.subplots(1, 2, figsize=(20, 10))
        # gray_next_obs = self.rgbd_to_grayscale(self.transformed_obs[0])[0, 0, :, :].numpy()
        # axes[0].imshow(gray_next_obs, cmap='gray', vmin=0, vmax=1.0)
        # transformed_rgb = np.swapaxes(np.swapaxes(self.transformed_obs[0].numpy()[:3, :, :], 0, -1), 0, 1)
        # axes[1].imshow(transformed_rgb)
        # # axes[1].matshow(mask)
        # plt.show()
        # plt.pause(0.001)
        return True
        # return self.ask_action_valid()

    def process_pixels(self, pix):
        post_pix = pix.copy().astype(np.float32)
        ratio_x = self.postcrop_pretransform_d.shape[0] / self.pretransform_depth.shape[0]
        ratio_y = self.postcrop_pretransform_d.shape[1] / self.pretransform_depth.shape[1]
        post_pix[0] *= ratio_x
        post_pix[1] *= ratio_y
        post_pix = post_pix.astype(np.uint16)
        post_pix[0] += WS_PC[0]
        post_pix[1] += WS_PC[2]

        return post_pix

    def check_pick_place_action(self, action_primitive, a_pick, a_place,
                                transformed_depth, transformed_rgb, rotation,
                                scale, value_map=None, pick_map=None,
                                all_value_maps=None, **kwargs):

        retval = super().check_pick_place_action(action_primitive, a_pick, a_place,
                                                 transformed_depth, transformed_rgb, rotation,
                                                 scale, value_map, pick_map, all_value_maps, **kwargs)

        valid_action = retval['valid_action']
        p1_grasp_cloth = retval['p1_grasp_cloth']

        if valid_action and p1_grasp_cloth:
            _ = self.visualize_and_ask_if_valid(retval)
            if self.number_of_restarts > 0:  # double check if the proposed action is valid
                valid_action = self.ask_action_valid()
        else:
            retval['valid_action'] = False
            # if self.number_of_restarts > 1:  # double check if the proposed action is valid
            #     _ = self.visualize_and_ask_if_valid(retval)
            return retval

        p1_cropped, p2_cropped = retval['pretransform_pixels'].copy()

        p1 = self.process_pixels(p1_cropped)
        p2 = self.process_pixels(p2_cropped)
        # Y seems to be a bit shifted
        # p1[0] -= 5  # 0 is X
        # p2[0] -= 5
        # p1[1] -= 5
        # p2[1] -= 5

        # transformed_pixels = np.array([a_pick[1:].numpy(), a_place[1:].numpy()])
        # pretransformed_pixels = np.array([p1, p2])

        # img = visualize_action_w_place(action_primitive=action_primitive, transformed_pixels=transformed_pixels,
        #                                pretransform_pixels=pretransformed_pixels,
        #                                value_map=value_map, pick_map=pick_map, all_value_maps=all_value_maps,
        #                                pretransform_depth=self.raw_pretransform_depth,
        #                                pretransform_rgb=self.raw_pretransform_rgb, transformed_rgb=transformed_rgb,
        #                                rotation=0, scale=1.0)
        #
        # im_real_dim = Image.fromarray(img)
        # im_real_dim.show()

        cam_intr = self.intrinsics_matrix

        y, x = p1
        np.save('depth_im.npy', self.raw_pretransform_depth.copy())  # TODO Remove, save just for testing
        np.save('raw_rgb.npy', self.raw_pretransform_rgb)  # TODO Remove, save just for testing
        # p1_grasp_cloth = self.preaction_mask[y, x]
        pick_point = list(pix_to_3d_position(x=x, y=y,
                                             depth_image=self.raw_pretransform_depth.copy(),
                                             cam_intr=cam_intr, cam_extr=self.top_cam_left_ur5_pose,
                                             cam_depth_scale=self.cam_depth_scale))
        y, x = p2
        p2_grasp_cloth = self.preaction_mask[y, x]
        place_point = list(pix_to_3d_position(x=x, y=y,
                                              depth_image=self.raw_pretransform_depth.copy(),
                                              cam_intr=cam_intr, cam_extr=self.top_cam_left_ur5_pose,
                                              cam_depth_scale=self.cam_depth_scale))

        retval.update({
            'valid_action': valid_action,
            'p1': pick_point,
            'p2': place_point,
            'p1_grasp_cloth': valid_action,
            'p2_grasp_cloth': p2_grasp_cloth
        })

        return retval

    def check_action(self, action_primitive, pixels,
                     transformed_depth, transformed_rgb,
                     scale, rotation,
                     value_map=None, all_value_maps=None,
                     **kwargs):
        # a_pick, a_place = pixels
        # kwargs.update({'a_pick': a_pick})
        # kwargs.update({'a_place': a_place})
        # retval = super().check_pick_place_action(**kwargs)
        retval = super().check_action(action_primitive, pixels,
                                      transformed_depth, transformed_rgb,
                                      scale, rotation,
                                      value_map=value_map, all_value_maps=all_value_maps, **kwargs)


        valid_action = retval['valid_action']
        p1_grasp_cloth = retval['p1_grasp_cloth']
        if valid_action and p1_grasp_cloth:
            _ = self.visualize_and_ask_if_valid(retval)
            if self.number_of_restarts > 0:  # double check if the proposed action is valid
                valid_action = self.ask_action_valid()
        else:
            retval['valid_action'] = False
            return retval

        p1_cropped, p2_cropped = retval['pretransform_pixels'].copy()

        p1 = self.process_pixels(p1_cropped)
        p2 = self.process_pixels(p2_cropped)
        # p1[0] -= 5  # 0 is X
        # p2[0] -= 5
        # p1[1] -= 5
        # p2[1] -= 5

        cam_intr = self.intrinsics_matrix

        y, x = p1
        # np.save('depth_im.npy', self.raw_pretransform_depth.copy())  # TODO Remove, save just for testing
        # np.save('raw_rgb.npy', self.raw_pretransform_rgb)  # TODO Remove, save just for testing
        # p1_grasp_cloth = self.preaction_mask[y, x]
        pick_point = list(pix_to_3d_position(x=x, y=y,
                                             depth_image=self.raw_pretransform_depth.copy(),
                                             cam_intr=cam_intr, cam_extr=self.top_cam_left_ur5_pose,
                                             cam_depth_scale=self.cam_depth_scale))
        y, x = p2
        p2_grasp_cloth = self.preaction_mask[y, x]
        place_point = list(pix_to_3d_position(x=x, y=y,
                                              depth_image=self.raw_pretransform_depth.copy(),
                                              cam_intr=cam_intr, cam_extr=self.top_cam_left_ur5_pose,
                                              cam_depth_scale=self.cam_depth_scale))

        retval.update({
            'valid_action': valid_action,
            'p1': pick_point,
            'p2': place_point,
            'p1_grasp_cloth': valid_action,
            'p2_grasp_cloth': p2_grasp_cloth
        })
        return retval

    def start_recording(self):
        if self.recording:
            return
        self.recording = True
        self.recording_daemon = setup_thread(target=self.record_video_daemon_fn)
        # self.recording_daemon_side = setup_thread(target=self.record_side_video_daemon_fn)

    def stop_recording(self):
        if not self.recording:
            return
        self.recording = False
        self.recording_daemon.join()
        # self.recording_daemon_side.join()
        self.recording_daemon = None
        # self.recording_daemon_side = None

    def record_side_video_daemon_fn(self):
        while self.recording:
            text = f'step {self.current_timestep}'
            if self.setup_cam is not None:
                if 'setup' not in self.env_video_frames:
                    self.env_video_frames['setup'] = []
                self.env_video_frames['setup'].append(
                    get_video_crop(self.setup_cam.get_rgbd(repeats=1)[0]))

    def record_video_daemon_fn(self):
        while self.recording:
            # NOTE: negative current_timestep's
            # act as error codes
            text = f'step {self.current_timestep}'
            if self.current_timestep == -1:
                text = 'randomizing cloth'
            elif self.current_timestep == -2:
                text = 'grasp failed'
            elif self.current_timestep == -3:
                text = 'cloth stuck'
            elif self.current_timestep == -4:
                text = 'ur5 timed out'
            if 'top' not in self.env_video_frames:
                self.env_video_frames['top'] = []
            top_view = cv2.resize(
                get_mask_workspace_crop(
                    self.top_cam.get_rgbd(repeats=1)[0].copy()),
                (256, 256))
            self.env_video_frames['top'].append(
                add_text_to_image(
                    image=top_view,
                    text=text, fontsize=16))
            if self.setup_cam is not None:
                if 'setup' not in self.env_video_frames:
                    self.env_video_frames['setup'] = []
                self.env_video_frames['setup'].append(
                    get_video_crop(self.setup_cam.get_rgbd(repeats=1)[0]))
            # if len(self.env_video_frames['setup']) > 50000:
            #     # episode typically ends in 4000 frames
            #     print('Robot probably got into error... Terminating')
            #     exit()

    def compute_coverage(self):
        coverage = compute_coverage(rgb=self.top_cam.get_rgbd()[0])
        # self.current_task.name[:-19], -19 as these are the date characters
        print(f"\tCoverage: {coverage/CLOTHS_DATASET[self.current_task.name[:-19]]['flatten_area']:.04f}")
        return coverage

    def get_obs(self):
        self.raw_pretransform_rgb, self.raw_pretransform_depth = self.top_cam.get_rgbd()

        # np.save(f"{self.root_dir}/raw_rgb_{self.current_timestep}.npy",
        #         self.raw_pretransform_rgb)

        self.postcrop_pretransform_rgb = get_workspace_crop(self.raw_pretransform_rgb.copy())
        self.postcrop_pretransform_d = get_workspace_crop(self.raw_pretransform_depth.copy())
        w, h = self.postcrop_pretransform_d.shape
        # assert w == h

        self.pretransform_rgb = cv2.resize(self.postcrop_pretransform_rgb, (self.image_dim, self.image_dim))
        self.pretransform_depth = cv2.resize(self.postcrop_pretransform_d, (self.image_dim, self.image_dim))

        # TODO Rotate image

        cloth_mask = self.get_cloth_mask(self.pretransform_rgb)


        # if not 'sequential' in self.action_primitive[0]:
        # if not self.action_primitive[0] == 'sequential-place' and not self.action_primitive[0] == 'sequential-dynamic':
        if self.replace_background:
            cloth_mask = (1 - cloth_mask).astype(bool)
            self.pretransform_rgb[..., 0][cloth_mask] = 0
            self.pretransform_rgb[..., 1][cloth_mask] = 0
            self.pretransform_rgb[..., 2][cloth_mask] = 0
        # if self.replace_background:
        #     # if self.action_primitive[0] == 'sequential-place':
        #     #     cloth_mask = cloth_mask.astype(bool)  # This works better for place
        #     # else:
        #         cloth_mask = (1-cloth_mask).astype(bool)
        #     self.pretransform_rgb[..., 0][cloth_mask] = 0
        #     self.pretransform_rgb[..., 1][cloth_mask] = 0
        #     self.pretransform_rgb[..., 2][cloth_mask] = 0

        # figure, axes = plt.subplots(1, 2, figsize=(20, 10))
        # axes[0].imshow(cloth_mask, cmap='gray', vmin=0, vmax=1.)
        # axes[1].imshow(self.pretransform_rgb)
        # #
        # plt.show()

        x, y = np.where(cloth_mask == 1)
        dimx, dimy = self.pretransform_depth.shape
        minx = x.min()
        maxx = x.max()
        miny = y.min()
        maxy = y.max()

        num_components = self._get_components_in_masked_rgb(self.pretransform_rgb)

        self.adaptive_scale_factors = self.scale_factors.copy()
        if self.use_adaptive_scaling:
            if self.compute_coverage()/CLOTHS_DATASET[self.current_task.name[:-19]]['flatten_area'] < 0.3:
                self.adaptive_scale_factors = self.adaptive_scale_factors[:4]

        if self.use_adaptive_scaling:
            try:
                # Minimum square crop
                cropx = max(dimx - 2*minx, dimx - 2*(dimx-maxx))
                cropy = max(dimy - 2*miny, dimy - 2*(dimy-maxy))
                crop = max(cropx, cropy)
                # Some breathing room
                crop = int(crop*1.5)
                if crop < dimx:
                    self.adaptive_scale_factors *= crop/dimx
                    self.episode_memory.add_value(
                        key='adaptive_scale',
                        value=float(crop/dimx))
            except Exception as e:
                print(e)
                print(self.current_task)
                exit()

        return preprocess_obs(self.pretransform_rgb.copy(), self.pretransform_depth.copy()), num_components

    def reset(self):
        self.episode_memory = SequentialMemory()
        self.episode_reward_sum = 0.
        self.current_timestep = -1
        self.terminate = False
        self.env_video_frames = {}
        # if self.dump_visualizations:        #     self.start_recording()
        self.current_task = self.get_task_fn()
        print("Moving out of the way")
        self.franka.out_of_the_way()
        print(f"New task={self.current_task.name}, difficulty = {self.current_task.task_difficulty}")
        if self.dump_visualizations:
            self.stop_recording()

        self.current_timestep = 0
        self.init_coverage = self.compute_coverage()

        obs, _ = self.get_obs()  # Check number of components or just
        self.episode_memory.add_value(key='pretransform_observations', value=obs)
        self.episode_memory.add_value(key='failed_grasp', value=0)
        self.episode_memory.add_value(key='timed_out', value=0)
        self.episode_memory.add_value(key='cloth_stuck', value=0)

        if self.randomise_rotation:
            random_rot_idx = np.random.randint(self.random_rotation.shape[0])
            self.rotations = [self.random_rotation[random_rot_idx]]

        self.transformed_obs = prepare_image(obs, self.get_transformations(), self.obs_dim,
                                             parallelize=self.parallelize_prepare_image)

        # Check if our policy uses grayscale or not
        if self.grayscale_only:
            gray_next_obs = self.rgbd_to_grayscale(self.transformed_obs[0])
            return gray_next_obs, self.ray_handle
        else:
            return self.transformed_obs, self.ray_handle

    def on_episode_end(self, log=False):
        if self.dump_visualizations and len(self.episode_memory) > 0:
            while True:
                hashstring = hashlib.sha1()
                hashstring.update(str(time()).encode('utf-8'))
                vis_dir = self.log_dir + '/' + hashstring.hexdigest()[:10]
                if not os.path.exists(vis_dir):
                    break
            os.mkdir(vis_dir)
            for key, frames in self.env_video_frames.items():
                if len(frames) == 0:
                    continue
                path = f'{vis_dir}/{key}.mp4'
                with imageio.get_writer(path, mode='I', fps=24) as writer:
                    for frame in (frames if not log else tqdm(frames, desc=f'Dumping {key} frames')):
                        writer.append_data(frame)
            # self.episode_memory.add_value(key='visualization_dir', value=vis_dir)
        self.env_video_frames.clear()

        self.stop_recording()
        super().on_episode_end(log=True)
        if os.path.exists(self.replay_buffer_path):
            with FileLock(self.replay_buffer_path + ".lock"):
                with h5py.File(self.replay_buffer_path, 'r') as file:
                    print('\tReplay Buffer Size:', len(file))
        print('='*10 + f'EPISODE END IN {self.current_timestep} STEPS' + '='*10)

    def check_action_reachability(self, **kwargs):
        return True, None

    def get_cam_pose(self):
        return self.top_cam_left_ur5_pose

    # PRIMITIVES TO MODIFY
    def pick_and_drag_primitive(self, p1, p2, move_speed: float, **kwargs):
        start_drag_pos, end_drag_pos = p1, p2
        start_drag_pos = bound_grasp_pos(start_drag_pos)
        end_drag_pos = bound_grasp_pos(end_drag_pos)

        if self.dump_visualizations:
            self.start_recording()

        predrag_point = deepcopy(start_drag_pos)
        end_predrag_point = deepcopy(end_drag_pos)
        print(f"Moving at speed={move_speed} Pick point ={start_drag_pos}, place point = {start_drag_pos}")

        predrag_point[2] += 0.2  # Z position - for going up need to
        end_predrag_point[2] += 0.2

        self.franka.send_home()

        self.franka.movel(params=[predrag_point + DEFAULT_ORN], traj_duration=5.0)

        self.franka.close_grippers_middle()

        # valid = self.ask_action_valid()
        # if not valid:
        #     self.franka.close_grippers_middle()
        #     self.franka.out_of_the_way()
        #     return valid

        self.franka.movel(params=[start_drag_pos + DEFAULT_ORN], traj_duration=2.0)
        self.franka.close_grippers()

        valid = self.ask_action_valid()
        if not valid:
            self.franka.close_grippers_middle()
            self.franka.movel(params=[predrag_point + DEFAULT_ORN], traj_duration=5.0)
            self.franka.out_of_the_way()
            if self.dump_visualizations:
                self.stop_recording()
            return valid

        self.franka.movel(params=[end_drag_pos + DEFAULT_ORN], traj_duration=move_speed)
        self.franka.close_grippers_middle()

        self.franka.movel(params=[end_predrag_point + DEFAULT_ORN], traj_duration=2.0)
        valid = self.ask_action_valid()
        if not valid:
            self.franka.out_of_the_way()
            if self.dump_visualizations:
                self.stop_recording()
            return valid
        # self.franka.open_grippers()

        self.franka.out_of_the_way()

        if self.dump_visualizations:
            self.stop_recording()

        return True

    def pick_and_place_primitive(self, p1, p2, move_speed: float, lift_height=0.2, **kwargs):
        pick_point, place_point = p1, p2
        pick_point = bound_grasp_pos(pick_point)
        place_point = bound_grasp_pos(place_point)
        if move_speed < 1.0:
            move_speed = 10.0

        if self.dump_visualizations:
            self.start_recording()

        print(f"Moving at speed={move_speed}, height={lift_height} Pick point ={pick_point}, place point = {place_point}")

        prepick_point = deepcopy(pick_point)
        prepick_point[2] += lift_height

        preplace_point = deepcopy(place_point)
        preplace_point[2] += lift_height

        self.franka.send_home()

        # valid = self.ask_action_valid()

        # if not valid:
        #     self.franka.out_of_the_way()
        #     raise GraspFailException

        # Move to the pre-pick point

        self.franka.movel(params=[prepick_point + DEFAULT_ORN], traj_duration=5.0)
        # self.franka.open_grippers()
        self.franka.close_grippers_middle()

        # valid = self.ask_action_valid()
        # if not valid:
        #     self.franka.close_grippers_middle()
        #     self.franka.out_of_the_way()
        #     return valid
        # Move to the pick position
        self.franka.movel(params=[pick_point + DEFAULT_ORN], traj_duration=2.5)

        self.franka.close_grippers()

        # valid = self.ask_action_valid()
        # if not valid:
        #     self.franka.close_grippers_middle()
        #     self.franka.out_of_the_way()
        #     return valid

        # valid = self.ask_action_valid()
        # if not valid:
        #     self.franka.open_grippers()
        #     self.franka.out_of_the_way()
        #     raise GraspFailException

        # self.ur5.movel(params=[backup_point + DEFAULT_ORN], j_vel=0.01, j_acc=0.01, blocking=True, use_pos=True)
        self.franka.movel(params=[prepick_point + DEFAULT_ORN], traj_duration=4.0)

        valid = self.ask_action_valid()
        if not valid:
            self.franka.close_grippers_middle()
            self.franka.out_of_the_way()
            if self.dump_visualizations:
                self.stop_recording()
            return valid

        self.franka.movel(params=[preplace_point + DEFAULT_ORN], traj_duration=move_speed)
        self.franka.movel(params=[place_point + DEFAULT_ORN], traj_duration=2.5)

        self.franka.close_grippers_middle()
        self.franka.movel(params=[preplace_point + DEFAULT_ORN], traj_duration=2.0)
        # self.franka.open_grippers()

        # if should_grasp_cloth and self.compute_iou() > 0.75:
        #     raise GraspFailException

        self.franka.out_of_the_way()
        if self.dump_visualizations:
            self.stop_recording()

        return True

    def dynamic_primitive(self, p1, p2, p1_grasp_cloth: bool, p2_grasp_cloth: bool,
                          mid_theta_ref: float, lift_height=0.2):
        start_dyn_pos, end_dyn_pos = p1, p2
        start_dyn_pos = bound_grasp_pos(start_dyn_pos)
        end_dyn_pos = bound_grasp_pos(end_dyn_pos)

        print(f"Moving at speed={mid_theta_ref} Pick point ={start_dyn_pos}, place point = {end_dyn_pos}")

        quintic_traj, tf = self.franka.compute_dynamic(init_pos=[start_dyn_pos],
                                                       end_pos=[end_dyn_pos],
                                                       mid_theta_ref=mid_theta_ref)

        if self.dump_visualizations:
            self.start_recording()

        # should_grasp_cloth = info['should_grasp_cloth']

        predyn_point = deepcopy(start_dyn_pos)
        end_predyn_point = deepcopy(end_dyn_pos)
        print(f"predrag_point point={predyn_point}")
        predyn_point[2] += lift_height
        end_predyn_point[2] += lift_height

        self.franka.send_home()

        self.franka.movel(params=[predyn_point + DEFAULT_ORN], traj_duration=5.0)

        self.franka.close_grippers_middle()

        # valid = self.ask_action_valid()
        # if not valid:
        #     self.franka.close_grippers_middle()
        #     self.franka.out_of_the_way()
        #     return valid

        self.franka.movel(params=[start_dyn_pos + DEFAULT_ORN], traj_duration=2.0)
        self.franka.close_grippers()

        valid = self.ask_action_valid()
        if not valid:
            self.franka.close_grippers_middle()
            self.franka.movel(params=[predyn_point + DEFAULT_ORN], traj_duration=5.0)
            self.franka.out_of_the_way()
            if self.dump_visualizations:
                self.stop_recording()
            return valid
        # Precompute it before executing then send the trajectory
        # Execute the dynamic primitive
        self.franka.movedynamic(quintic_traj=quintic_traj, tf=tf)
        valid = self.ask_action_valid()
        if not valid:
            self.franka.close_grippers_middle()
            self.franka.movel(params=[predyn_point + DEFAULT_ORN], traj_duration=5.0)
            self.franka.out_of_the_way()
            if self.dump_visualizations:
                self.stop_recording()
            return valid
        self.franka.movel(params=[end_dyn_pos + DEFAULT_ORN], traj_duration=1.0)
        self.franka.close_grippers_middle()

        self.franka.movel(params=[end_predyn_point + DEFAULT_ORN], traj_duration=2.0)

        # self.franka.open_grippers()

        self.franka.out_of_the_way()

        if self.dump_visualizations:
            self.stop_recording()

        return True

    def rgbd_to_grayscale(self, rgbd_image):
        "Modify the method so that it also inverts the gray image"
        image_transposed = np.transpose(rgbd_image[:3, :, :].numpy(), (1, 2, 0))
        # RGB-> GRAY Shape (H x W) -> [1, 1, H, W]
        grayscale_img = cv2.cvtColor(image_transposed, cv2.COLOR_RGB2GRAY)

        # inverted_gray_img = np.abs(grayscale_img - 1)

        # figure, axes = plt.subplots(1, 2, figsize=(20, 10))
        # axes[0].imshow(grayscale_img, cmap='gray', vmin=0, vmax=1)
        # axes[0].set_title("Grayscale raw")
        # axes[1].imshow(inverted_gray_img, cmap='gray', vmin=0, vmax=1)
        # axes[1].set_title("Inverted grayscale")
        # # axes[1].matshow(mask)
        # plt.show()
        # plt.pause(0.001)

        end_img = grayscale_img[None, None, :, :]
        # end_img = inverted_gray_img[None, None, :, :]

        # Could modify the background and make it a bit gray
        # end_img[]

        return torch.from_numpy(end_img).to(rgbd_image.device)
        # return torch.from_numpy(grayscale_img).to(rgbd_image.device)

