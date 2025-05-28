resFilename = 'rs.txt';
gtFilename = 'gt.txt';

% Read gt file
gtdata = dlmread(gtFilename);

% Read result file
resdata = dlmread(resFilename);

% Evaluate sequence
threshold = 2.73e-3;
[recall, precision, idswitches, MOTA, MOTP] = CLEAR_MOT(gtdata, resdata, threshold);

fprintf('MOTA: %.2f\nMOTP: %.2f\n',  MOTA, MOTP);