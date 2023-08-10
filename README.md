# dynamics_bag
- In addition to David’s comments: change to remote_franka3 in both header files of Franka 3

- Enable commments in .bashrc when running dual franka to get command from franka2


For runnning single robot:
- rosrun franka_david frankapy


For running dual robots:
- from franka2: roscore, "rosrun franka_david remotefranka"
- enable lines in .bashrc in franka3, then "rosrun franka_david remotefranka"
- then "python script eric_test_dual_franka.py" in franka2
- rosrun franka_david remotefranka

- git add ".gitignore" and push it after editing ignore rules, but before they take effect

For OpiTrack capture:
- put two markers close enough and they can be seen as a single point visible from both the inside and outisde of the rim (helps with occlusions)
    - usually seen as a single point in OptiTrack, but do additional filtering to combine monitored points that are captured too close together so that they will not
    -   cause convex hull calculation to fail!
