function disp_sim_results(L, R, d, GT)
% function for displaying simulation results
% disp_sim_results(L, R, d, GT)
%
% L, R are the left and right image pair
% d is the estimated disparity map
% GT is the ground truth

% display images side by side
stdisp(L, R);

% display bar graph side-by-side comparisson of disparities to ground truth 
figure
idisp([d, GT], 'bar');

end