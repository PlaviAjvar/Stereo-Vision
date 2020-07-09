function [L, R, d, GT, sim, DSI] = sim_method(image, method)
% Utility function for simulating generic stereo algorithm 
% [L, R, d, GT, sim, DSI] = sim_method(image, method)
%
% image = {'Aloe', 'Lampshade1', 'Rocks1'}
% method = {'Baseline', 'Classic', 'SmoothDP', 'OrderDP', 'SGM', 'LoopyBP', 'MultiscaleDP'}
%
% L, R are the image pair
% d is the disparity map
% GT is the ground truth
% sim is the similarity measure for matching pixel with given disparity
% DSI is the disparity space image


% load images and display
L = iread(strcat(image, '\view1.png'), 'reduce', 10);
R = iread(strcat(image, '\view5.png'), 'reduce', 10);

% load ground truth and adjust for size
GT = iread(strcat(image,'\disp1.png'), 'reduce', 10);
GT = GT / 2;
dmin = double(min(min(GT)));
dmax = double(max(max(GT)));

% fix correlation window and call utility function for finding disparities
H = 5;
[d, sim, DSI] = find_disparity(L, R, [dmin, dmax], H, method);

% d = istereo(L, R, [dmin, dmax], H);

% load image offset (cropped pixels from left image) and adjust disparities
fileID = fopen(strcat(image,'\dmin.txt'));
offset = fscanf(fileID, '%d') / 2;
d = d + offset;
GT = GT + offset;

end