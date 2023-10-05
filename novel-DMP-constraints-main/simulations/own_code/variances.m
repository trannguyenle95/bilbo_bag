clear all
clc

W1 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d1_flip.csv'));
W2 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d2_flip.csv'));
W3 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d3_flip.csv'));
W4 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d4_flip.csv'));
W5 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d5_flip.csv'));
W6 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d6_flip.csv'));
W7 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d7_flip.csv'));
W8 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d8_flip.csv'));
W9 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d9_flip.csv'));
W10 = csvread(fullfile('/home/erichannus/catkin_ws/src/Data/Weights/W_d10_flip.csv'));


%WIP - need separate analysis for each DoF? 
% OR flatten each demo to get 7*50 long vector and then calculate cov
% matix for this?? << TRY THIS
%Cov_j1 = W1(1,:)
%V0 = var(cat(3,W1,W2, W3),0,3);  
%V1 = var(cat(3,W1,W2, W3),1,3);  

W_variances = var(cat(3,W1,W2, W3),0,3);  
imagesc(W_variances)
colorbar