# Handbook

### 1. Prepare
### step1 - 张正友相机标定
toolkit: 张正友
input: calib1, calib2
output: calib.mat

### step2 - 标定数据格式转换
toolkit: formatKit/calib2stereModel.m
input: calib.mat
output: stereoModel.mat

### step3 - 数据格式转换
toolkit: formatKit/video2png.m
input: 1.wmv, 2.wmv
output: cam1/fr%05.png, cam2/fr%05.png

### step4 - 盒子计算
toolkit: formatKit/calculateCube.m
input: 视图1,2上手动标的顶点, stereoModel.mat
output: cube.mat

## 2 detect & track

### step5 - 检测跟踪
toolkit: 主程序
input: cam1/fr%05.png, cam2/fr%05.png
cube.mat, stereoModel.mat, flyModel.mat
output: result.mat (主要是tk)

## 3 post-progress

### step6 - 可视化
toolkit: analyseKit
input: result.mat
output: 

