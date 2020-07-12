function good_percentage = postprocess(L, R, d, sim, DSI, GT)
% Function for postprocessing an image, displaying relevant information
% Works only for Classic approach
% good_percentage = postprocess(L, R, d, DSI)
%
% L, R are the image pair
% d is the disparity image
% DSI is the disparity space image
% sim is the similarity score
% GT is the ground truth
%
% good_percentage is the percentage of pixels which pass the tests

% taken from Corke
status = zeros(size(d));
[U,V] = imeshgrid(L);

% remove NaN values
status(isnan(d)) = 3;

% (C) ignore peek if it's energy is too high (low similarity)
e_high = 0.4;
status(sim > e_high) = 1;

% (D) interpolate polynomial through points next to optimum
% check quadratic coefficient
flat = 0.05;
[A, B, C] = quad_int(DSI);
status(A < flat) = 2;

% display on image
idisp(status);
colormap(colorname({'lightgreen', 'red', 'blue', 'orange'}));

% percentage of "green" samples (good matches)
good_percentage = sum(status(:) == 0) / prod(size(status)) * 100;

% invalidate bad matches in disparity map
d(status > 0) = NaN;

% display unreliable results
figure
ipixswitch(isnan(d), 'red', d/90);

end