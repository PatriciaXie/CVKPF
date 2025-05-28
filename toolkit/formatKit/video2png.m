clear; clc; close all;

nFly = 1; Case = 1;

video1 = 'D:\02-dataset\xie2019\3D\20200109\exp1\exp1_C1.wmv'; 
video2 = 'D:\02-dataset\xie2019\3D\20200109\exp1\exp1_C2.wmv'; 

path = sprintf('data/realdata%03d-%d', nFly, Case);
path1 = sprintf('data/realdata%03d-%d/cam1', nFly, Case);
path2 = sprintf('data/realdata%03d-%d/cam2', nFly, Case);
mkdir(path);
mkdir(path1);
mkdir(path2);

v1 = VideoReader(video1);
v2 = VideoReader(video2);
k = 1; T = 5000; % T = v1.Duration * v1.FrameRate;
waitHandle=waitbar(0, 'ÊÓÆµ×ª»»µ½png...');
while k <= T
    img1 = rgb2gray(readFrame(v1)); 
    img2 = rgb2gray(readFrame(v2)); 
    filename1 = sprintf('%s/fr%05d.png', path1, k);
    filename2 = sprintf('%s/fr%05d.png', path2, k);
    imwrite(img1, filename1);
    imwrite(img2, filename2);
    waitbar(k/T);
     k = k + 1;
end
close(v1);
close(v2);
delete(waitHandle);
clear waitHandle