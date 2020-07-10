function dif = dif_metric(d, GT)
% Function for calculating mean absolute error between images
% dif = dif_metric(d, GT)
%
% d is the disparity image (the algorithm output)
% GT is the ground truth

delta = double(d) - double(GT);
entries = ~isnan(delta);
dif = sum(abs(delta(entries))) / sum(entries(:));

end