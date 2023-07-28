import sys
sys.path.append("/home/david/catkin_ws/src/franka_david/scripts/qdc-manip")

import time
import numpy as np
import cv2
from scipy import optimize
from copy import deepcopy

from kinect import KinectClient
from franka import Franka, DEFAULT_ORN, bound_grasp_pos


def pick_and_place_primitive(franka, p1, p2, height=0.2, **kwargs):
    pick_point, place_point = p1, p2
    pick_point = bound_grasp_pos(pick_point)
    place_point = bound_grasp_pos(place_point)

    print(f"Pick point ={pick_point}, place point = {place_point}")

    prepick_point = deepcopy(pick_point)
    backup_point = deepcopy(pick_point)
    # - for going up
    prepick_point[2] += 0.1
    preplace_point = deepcopy(place_point)
    preplace_point[2] += 0.1

    # valid = self.ask_action_valid()

    # if not valid:
    #     self.ur5.out_of_the_way()
    #     raise GraspFailException

    # Move the franka to the center of the table with joints positions
    franka.send_home()

    # Move to the pre-pick point
    franka.movel(params=[prepick_point + DEFAULT_ORN], traj_duration=5.0)
    franka.open_grippers()
    # valid = self.ask_action_valid()
    #
    # if not valid:
    #     self.ur5.out_of_the_way()
    #     raise GraspFailException

    franka.open_grippers_middle()

    # Move to the pick position
    franka.movel(params=[pick_point + DEFAULT_ORN], traj_duration=1.0)

    franka.close_grippers()

    # valid = self.ask_action_valid()
    # if not valid:
    #     self.ur5.open_grippers()
    #     self.ur5.out_of_the_way()
    #     raise GraspFailException

    # self.ur5.movel(params=[backup_point + DEFAULT_ORN], j_vel=0.01, j_acc=0.01, blocking=True, use_pos=True)
    franka.movel(params=[prepick_point + DEFAULT_ORN])
    # valid = self.ask_action_valid()
    # if not valid:
    #     self.ur5.open_grippers()
    #     self.ur5.out_of_the_way()
    #     raise GraspFailException
    franka.movel(params=[preplace_point + DEFAULT_ORN])
    franka.movel(params=[place_point + DEFAULT_ORN], traj_duration=1.0)
    franka.open_grippers_middle()
    franka.movel(params=[preplace_point + DEFAULT_ORN], traj_duration=2.0)
    franka.open_grippers()


if __name__ == '__main__':

    z_position = 0.125
    bottom_left = [[0.30053, -0.35026, z_position]]
    bottom_right = [[0.72368, -0.35026, z_position]]  # 78265126
    top_left = [[0.30053, 0.35026, z_position]]
    top_right = [[0.72368, 0.35026, z_position]]  # 7761067
    # pick_point = [0.38270446324518415, -0.3335766877080274, 0.11399022352134369]
    # place_point = [0.30650399048670585, -0.2497963299321307, 0.11918521773636348]
    pick_point = [0.4641766529428322, 0.3429609653479405, 0.15]
    place_point = [0.32, -0.32, 0.15]
    # pick_point = [0.47874473539475526, -0.32809403611756754, 0.11379303542523656]
    # place_point = [0.32988289641969876, -0.2556129424156895, 0.12097379715571033]

    robot = Franka(init_node=True)
    # input('Close grippers')
    # robot.close_grippers()
    # input('Open grippers')
    # robot.open_grippers()
    # robot.out_of_the_way()

    # p1 = [[bottom_left[0][0] + 0.1, bottom_left[0][1] + 0.1, 0.15]]  # 4.6s
    # p2 = [[bottom_left[0][0] + 0.05, bottom_left[0][1] + 0.05, 0.15]]  # 2.3s
    # p3 = [[bottom_left[0][0] + 0.02, bottom_left[0][1] + 0.02, 0.15]]  # 0.9s
    # p4 = [[bottom_left[0][0] + 0.02, bottom_left[0][1], 0.15]]

    # robot.send_home()
    # robot.movel(bottom_left, traj_duration=4.0)
    # input('Next Movement')
    # robot.movel(bottom_right, traj_duration=4.0)
    # input('Next Movement')
    # robot.send_home()
    # robot.movel(top_left, traj_duration=4.0)
    # input('Next Movement')
    # robot.movel(top_right, traj_duration=4.0)
    # input('Next Movement')

    # robot.send_home()
    # robot.movel(bottom_left, traj_duration=4.0)
    quintic_traj, tf = robot.compute_dynamic(init_pos=bottom_left, end_pos=top_right,
                                             mid_theta_ref=-0.5)
    # max velocity -0.5 - min velocity -0.08
    # -0.05 -- 34 seconds
    # -0.08 -- 26 seconds
    # mid_theta_ref = -0.1 -- 20seconds
    # -0.2 -- 11 seconds  # Not that fast actually
    # 30 seconds?
    # input("Perform dynamic movement")
    # robot.movedynamic(quintic_traj=quintic_traj, tf=tf)



    # input("Test next trajectory")
    #
    # robot.send_home()
    # robot.movel(bottom_left, traj_duration=4.0)
    # robot.movedynamic(bottom_left, p2, mid_theta_ref=-0.1)
    # input("Test next trajectory")
    #
    #
    # robot.send_home()
    # robot.movel(bottom_left, traj_duration=4.0)
    # robot.movedynamic(bottom_left, p3, mid_theta_ref=-0.1)
    # input("Test next trajectory")
    #
    # robot.send_home()
    # robot.movel(bottom_left, traj_duration=4.0)
    # robot.movedynamic(bottom_left, p4, mid_theta_ref=-0.1)
    # input("Test next trajectory")

    # robot.send_home()

    # From bottom left to bottom right
    # mid theta ref | tf (s) -- K = 0.005, does not find solution, final tf found
    # -1.2 | 0.36
    # -0.1 | 3.83
    # -0.05 | 5.15
    # -0.01 | 6.39
    # -0.005 | 8.12
    # mid theta ref | tf (s) -- K = 0.05, does not find solution, final tf found
    # -0.01 | 20 seconds -- might be quasi-static enough
    # -0.1 | 20 seconds / 17 seconds


    a = 0
