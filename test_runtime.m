function runtime = test_runtime(image, method, H)
% Function which tests runtime by repeatedly running algorithm and taking
% average. Default number of times run is 5.
% runtime = test_runtime(image, method, H)
%
% image is the subfolder containing Middlebury image
% image = {'Aloe', 'Lampshade1', 'Rocks1'}
%
% method is a string describing the method in use
% method = {'Baseline', 'Classic', 'SmoothDP', 'OrderDP', 'SGM', 'LoopyBP'}
%
% H is the halfsize of the correlation window
% runtime is the execution time in miliseconds

% how many repeats
no_try = 5;
total_time = 0;

for iter = 1:no_try
    [L, R, d, GT, sim, DSI, exec_time] = sim_method(image, method, H);
    total_time = total_time + exec_time;
end

runtime = total_time / no_try;

end