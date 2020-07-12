# Stereo-Vision

Implementing various dense stereo matching algorithms. Most of the algorithms are implemented in C, because of speed. The image processing is handled by MATLAB.

## Install
The requirements are Peter Corke's machine vision toolbox and the MinGW-w64 C compiler. Install instructions for the former can be found on Corke's website. The second can be installed directly through MATLAB.

## Functionality
Documentation for the functions can also be found through MATLAB's help.

>  function [L, R, d, GT, sim, DSI] = sim_method(image, method)

Simulates function given by parameter "method", for given image subfolder "image". The possible values for "image" are {'Aloe', 'Lampshade1', 'Rocks1'}, although other folders can be added for testing. The method can be one of the following: {'Baseline', 'Classic', 'SmoothDP', 'OrderDP', 'SGM', 'LoopyBP'}.

Baseline calls Corke's matching algorithm, based on ZNCC. Classic uses the same approach, but my implementation. SmoothDP implements scanline dynamic programming, with smoothness costs. Order DP implements ordering constrained DP. SGM implements semi-global matching. LoopyBP implements serial loopy belief propagation.

The output parameters are the image pair "L, R" (half-size version), the disparity image "d", the ground truth "GT", the similarity measure map "sim", and the disparity space image "DSI" (3D representation).

> function d_new = quadratic_process(d, DSI, offset)

Uses quadratic interpolation to obtain better estimates for disparities. Offset represents the value to be added to the disparities, due to image cropping.

> function [A, B, C] = quad_int(DSI)

Implements quadratic interpolation. "A, B, C" are matrices of values, representing coefficients for each pixel.

> function good_percentage = postprocess(L, R, d, sim, DSI, GT)

Analyzes percentage of good matches (good_percentage). Also displays image which displays problematic pixels.

> function [Mlam, Msat, Mmae] = param_sim(image, method)

Used for parameter simulation. Specify image and method, and in the table Mmae and Mmse contain the errors for various parameter choices.

> function dif = mse_metric(d, GT)
> function dif = mae_metric(d, GT)

Implement the error metrics between disparity image and ground truth.

> function [d, sim, DSI, exec_time] = find_disparity(L, R, d_range, H, method)

Call optimization methods directly, rather than through sim_method. Very similar in structure. The exec_time is the execution time of the algorithm in miliseconds.

> function disp_results(L, R, d, GT)

Displays the image pair L, R side by side using Corke's stdisp. Also displays disparity image and ground truth on a bar graph, side by side.