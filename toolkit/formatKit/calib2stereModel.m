clear; clc; close all;
% load('D:\02-dataset\chen2017\stStereoModel.mat');
calib = load('data\realdata001-1\calib.mat'); 
cam1 = calib.stereoParams.CameraParameters1;
cam2 = calib.stereoParams.CameraParameters2;
stereo = calib.stereoParams;

stereoModel.camResolution = cam1.ImageSize;
stereoModel.fundamentals{1} = stereo.FundamentalMatrix;
stereoModel.cams(1).projection = cam1.IntrinsicMatrix' * [eye(3), zeros(3, 1)];
stereoModel.cams(1).intrinsic = cam1.IntrinsicMatrix';
stereoModel.cams(1).center = zeros(1, 3);

R = stereo.RotationOfCamera2';
T = stereo.TranslationOfCamera2;
stereoModel.cams(2).projection = cam2.IntrinsicMatrix' * [R, T']; % 以相机1为世界系
stereoModel.cams(2).intrinsic = cam2.IntrinsicMatrix';
stereoModel.cams(2).center = T;

save('data\realdata001-1\stereoModel.mat', 'stereoModel');