function [d, sim, DSI, exec_time] = find_disparity(L, R, d_range, H, method)
% Helper function for calling algorithm of choice
% [d, sim, DSI, exec_time] = find_disparity(L, R, d_range, H, method)
%
% L, R are the image pair
% d_range = [d_min, d_max] is the range of feasible disparities
% H is the half-length of the correlation window
% method = {'Baseline', 'Classic', 'SmoothDP', 'OrderDP', 'SGM', 'LoopyBP', 'MultiscaleDP'}

if strcmp(method, 'Baseline')
    [d, sim, DSI] = istereo(L, R, d_range, H);
else
    % ensure images are monochromatic
    L = imono(L);
    R = imono(R);
    
    % call C solver
    [DSI, exec_time] = fast_stereo(L, R, H, d_range, method);
    
    % get disparity estimates from disparity state image
    [sim, d] = min(DSI, [], 3);
    d(isnan(sim)) = NaN;

    % adjust for offset
    d = d + (d_range(1) - 1);
end

end