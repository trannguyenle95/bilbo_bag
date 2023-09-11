function robot = createFranka()
%https://se.mathworks.com/help/robotics/ref/rigidbodytree.html

%https://frankaemika.github.io/docs/control_parameters.html

%NOTE franka gives DH parameters in Craig's convention (modified),
%so use "mdh" argument in https://se.mathworks.com/help/robotics/ref/rigidbodyjoint.setfixedtransform.html
% and change order given in franka source to match [a alpha d theta] used
% in Matlab:


dhparams = [0           0           0.333        0;
            0	        -pi/2       0            0;
            0	        pi/2        0.316        0;
            0.0825      pi/2        0            0;
            -0.0825     -pi/2       0.384        0;
            0           pi/2        0            0;
            0.088       pi/2        0.107+0.1    0]; %add offset to flange plus 10cm offset to get to center of hand


%seems to be same DH params for Research 3: https://www.generationrobots.com/media/franka-emika-research-3-robot-datasheet.pdf

%last row "theta" is ignored, so can be set to zero for all

robot = rigidBodyTree;

body1 = rigidBody('body1');
jnt1 = rigidBodyJoint('jnt1','revolute');
body2 = rigidBody('body2');
jnt2 = rigidBodyJoint('jnt2','revolute');
body3 = rigidBody('body3');
jnt3 = rigidBodyJoint('jnt3','revolute');
body4 = rigidBody('body4');
jnt4 = rigidBodyJoint('jnt4','revolute');
body5 = rigidBody('body5');
jnt5 = rigidBodyJoint('jnt5','revolute');
body6 = rigidBody('body6');
jnt6 = rigidBodyJoint('jnt6','revolute');
body7 = rigidBody('body7');
jnt7 = rigidBodyJoint('jnt7','revolute');

%https://se.mathworks.com/help/robotics/ref/rigidbodyjoint.setfixedtransform.html
%NOTE changed to use modified DH params
%Given in order [a alpha d theta]
setFixedTransform(jnt1,dhparams(1,:),'mdh');
setFixedTransform(jnt2,dhparams(2,:),'mdh');
setFixedTransform(jnt3,dhparams(3,:),'mdh');
setFixedTransform(jnt4,dhparams(4,:),'mdh');
setFixedTransform(jnt5,dhparams(5,:),'mdh');
setFixedTransform(jnt6,dhparams(6,:),'mdh');
setFixedTransform(jnt7,dhparams(7,:),'mdh');

body1.Joint = jnt1;
body2.Joint = jnt2;
body3.Joint = jnt3;
body4.Joint = jnt4;
body5.Joint = jnt5;
body6.Joint = jnt6;
body7.Joint = jnt7;

%Add constraints to joints
%https://se.mathworks.com/help/robotics/ref/rigidbodyjoint.html:
%CURRENTLY BREAKS IK - why?
%HOWEVER, setting these limits makes show(franka) result look better
body1.Joint.PositionLimits = [-2.7437 2.7437];
body2.Joint.PositionLimits = [-1.7628 1.7628];
body3.Joint.PositionLimits = [-2.8973 2.8973];
body4.Joint.PositionLimits = [-0.1518 -3.0421];
body5.Joint.PositionLimits = [-2.8065 2.8065];
body6.Joint.PositionLimits = [0.5445 3.7525];
body7.Joint.PositionLimits = [-2.8973 2.8973];

addBody(robot,body1,'base');
addBody(robot,body2,'body1');
addBody(robot,body3,'body2');
addBody(robot,body4,'body3');
addBody(robot,body5,'body4');
addBody(robot,body6,'body5');
addBody(robot,body7,'body6');

robot.DataFormat = 'column'; %so vector can be used instead of struct in inverse kinematics

%{
robot = importrobot("frankaEmikaPanda.urdf")
removeBody(robot, 'panda_rightfinger')
removeBody(robot, 'panda_leftfinger')

robot.DataFormat = 'column'; %so vector can be used instead of struct in inverse kinematics
%}


end