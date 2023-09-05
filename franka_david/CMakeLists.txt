cmake_minimum_required(VERSION 3.4)
project(franka_david)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(catkin REQUIRED COMPONENTS
  controller_interface
  dynamic_reconfigure
  eigen_conversions
  franka_hw
  franka_gripper
  geometry_msgs
  hardware_interface
  joint_limits_interface
  tf
  tf_conversions
  message_generation
  pluginlib
  realtime_tools
  roscpp
  rospy
  urdf
  visualization_msgs
)

find_package(Eigen3 REQUIRED)
find_package(Franka 0.9.0 QUIET)
if(NOT Franka_FOUND)
  find_package(Franka 0.8.0 REQUIRED)
endif()

add_message_files(FILES
MotionPy.msg
MovetoPy.msg
GripperPy.msg
GripperGraspPy.msg
JointMotionPy.msg
JointTrajPy.msg
RemoteMotionPy.msg
)

generate_messages(DEPENDENCIES geometry_msgs)

catkin_package(
  INCLUDE_DIRS include cpptoolbox/include/
  LIBRARIES cpp_fares
  CATKIN_DEPENDS
    controller_interface
    dynamic_reconfigure
    eigen_conversions
    franka_hw
    franka_gripper
    geometry_msgs
    hardware_interface
    joint_limits_interface
    tf
    tf_conversions
    message_runtime
    pluginlib
    realtime_tools
    roscpp
    urdf
    visualization_msgs
  DEPENDS Franka
)

add_library(cpp_fares
cpptoolbox/include/Exceptions.h
cpptoolbox/include/Frame3D.h
cpptoolbox/include/include.h
cpptoolbox/include/mathUtility.hpp
cpptoolbox/include/matrix.hpp
cpptoolbox/include/numericalAlg.h
cpptoolbox/include/quaternion.hpp
cpptoolbox/include/utility_io.h
cpptoolbox/src/columnvector.cpp
cpptoolbox/src/Frame3D.cpp
cpptoolbox/src/main.cpp
cpptoolbox/src/mathUtility.cpp
cpptoolbox/src/matrix3d.cpp
cpptoolbox/src/matrix.cpp
cpptoolbox/src/matrixutils.cpp
cpptoolbox/src/numericalAlg.cpp
cpptoolbox/src/PositionVector3D.cpp
cpptoolbox/src/quaternion.cpp
cpptoolbox/src/RotationMatrix.cpp
cpptoolbox/src/rowvector.cpp
cpptoolbox/src/Transform3D.cpp
cpptoolbox/src/Twist.cpp
cpptoolbox/src/utility_io.cpp
cpptoolbox/src/Wrench.cpp
# Include the examples common as part of the library
include/examples_common.h
src/examples_common.cpp
)

add_dependencies(cpp_fares
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

target_link_libraries(cpp_fares PUBLIC
  ${Franka_LIBRARIES}
  ${catkin_LIBRARIES}
)

target_include_directories(cpp_fares SYSTEM PUBLIC
  ${Franka_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${catkin_INCLUDE_DIRS}
)

target_include_directories(cpp_fares PUBLIC
  include
  cpptoolbox/include
)

catkin_install_python(PROGRAMS
  #scripts/simple_ros.py 
  scripts/test_franka.py
  scripts/calibrate_camera_franka_arm.py
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)


add_executable(frankapy
  src/control_through_python.cpp
)
target_include_directories(frankapy PUBLIC
  ${catkin_INCLUDE_DIRS}
  cpptoolbox/include/
  include
)
target_link_libraries(frankapy PUBLIC
  ${Franka_LIBRARIES}
  ${catkin_LIBRARIES}
  cpp_fares
)
add_dependencies(frankapy
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

add_executable(remotefranka
  src/control_remote.cpp
)
target_include_directories(remotefranka PUBLIC
  ${catkin_INCLUDE_DIRS}
  cpptoolbox/include/
  include
)
target_link_libraries(remotefranka PUBLIC
  ${Franka_LIBRARIES}
  ${catkin_LIBRARIES}
  cpp_fares
)
add_dependencies(remotefranka
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

