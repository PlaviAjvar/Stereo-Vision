function [L, R, d, GT, sim, DSI, exec_time] = sim_method(image, method, H)
% Utility function for simulating stereo algorithms
% [L, R, d, GT, sim, DSI, exec_time] = sim_method(image, method, H)
%
% image = {'Aloe', 'Lampshade1', 'Rocks1'}
% method = {'Baseline', 'Classic', 'SmoothDP', 'OrderDP', 'SGM', 'LoopyBP'}
% H is the halfsize of the correlation window
%
% L, R are the image pair
% d is the disparity map
% GT is the ground truth
% sim is the similarity measure for matching pixel with given disparity
% DSI is the disparity space image
% exec_time is the execution time given in miliseconds


% load images
L = iread(strcat(image, '\view1.png'));
R = iread(strcat(image, '\view5.png'));

% load ground truth and adjust for size
GT = iread(strcat(image,'\disp1.png'));
GT = GT / 2; 
dmin = double(min(min(GT)));
dmax = double(max(max(GT)));

% fix correlation window and call utility function for finding disparities
[d, sim, DSI, exec_time] = find_disparity(L, R, [dmin, dmax], H, method);

% d = istereo(L, R, [dmin, dmax], H);

% load image offset (cropped pixels from left image) and adjust disparities
fileID = fopen(strcat(image,'\dmin.txt'));
offset = fscanf(fileID, '%d') / 2;
d = d + offset;
GT = GT + offset;

end