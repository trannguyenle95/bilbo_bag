import os
import rospy
from std_msgs.msg import Bool
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('CurrentFranka', type=int, help='Franka version of CURRENT robot: 2 or 3')
    args = parser.parse_args()    
    
    datafolder = os.path.join(os.path.expanduser('~'), 'catkin_ws', 'src', 'Data')
    if args.CurrentFranka == 2:
        filepath = os.path.join(datafolder+"/"+"fr3_error.txt") #read errors from other robot
    elif args.CurrentFranka == 3:
        filepath = os.path.join(datafolder+"/"+"fr2_error.txt") #read errors from other robot
    else:
        raise Exception("Franka must be either 2 or 3") 

    if os.path.exists(filepath):
        os.remove(filepath) 

    def callback(data):
        print("callback for stopping")
        if data.data == True:
            f = open(filepath, 'x')

    rospy.init_node('listener', anonymous=True)
    if args.CurrentFranka == 2:
        rospy.Subscriber("/franka3_control_node/franka2_state", Bool, callback) #read errors from other robot
    elif args.CurrentFranka == 3:
        rospy.Subscriber("/franka2_control_node/franka2_state", Bool, callback) #read errors from other robot

    rospy.spin()