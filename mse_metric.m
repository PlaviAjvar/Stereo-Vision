function dif = mse_metric(d, GT)
% Function for calculating mean square error between disparity image and GT
% dif = mae_metric(d, GT)
%
% d is the disparity image (the algorithm output)
% GT is the ground truth

delta = double(d) - double(GT);
entries = ~isnan(delta);
dif = sum(delta(entries) .* delta(entries)) / sum(entries(:));

end