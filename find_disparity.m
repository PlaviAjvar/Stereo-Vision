function [d, sim, DSI] = find_disparity(L, R, d_range, H, method)
% Helper function for calling algorithm of choice
% [d, sim, DSI] = find_disparity(L, R, d_range, H, method)

% L, R are the image pair
% d_range = [d_min, d_max] is the range of feasible disparities
% H is the half-length of the correlation window
% method = {'Baseline', 'Classic', 'SmoothDP', 'OrderDP', 'SGM', 'LoopyBP', 'MultiscaleDP'}

if strcmp(method, 'Baseline')
    [d, sim, DSI] = istereo(L, R, d_range, H);
end

end