function [Mlam, Msat, Mmae, Mmse] = param_sim(image, method)
% Function used to simulate parametrized algorithms
% [Mlam, Msat, Mmae, Mmse] = param_sim(image, method)
%
% method = {'ScanlineDP', 'SGM', 'LoopyBP'}
%
% Mlam, Msat are the matrices of scalar/saturation pairs (meshgrid)
% Mmae is the mean absolute error for GT and d
% Mmse is the mean square error for GT and d

lam  = [0.05 0.1 0.2 0.5];
vsat = [0.5 0.75 1  2];

[Mlam, Msat] = meshgrid(lam, vsat);
    
% load images
L = iread(strcat(image, '\view1.png'));
R = iread(strcat(image, '\view5.png'));
L = imono(L);
R = imono(R);

% load ground truth and adjust for size
GT = iread(strcat(image,'\disp1.png'));
GT = GT / 2; 
dmin = double(min(min(GT)));
dmax = double(max(max(GT)));
H = 5;
    
Mmae = zeros(size(lam,2));
Mmse = zeros(size(lam,2));
fileID = fopen(strcat(image,'\dmin.txt'));
offset = fscanf(fileID, '%d') / 2;
GT = GT + offset;
    
for i = 1:size(lam,2)
    for j = 1:size(vsat,2)
        if Mlam(i,j) < Msat(i,j)
            DSI = fast_stereo(L, R, H, [dmin, dmax], method, [Mlam(i,j), Msat(i,j)]);
            [sim, d] = min(DSI, [], 3);
            d(isnan(sim)) = NaN;
            d = d + offset;
            d = d + (dmin - 1);
            Mmae(i,j) = mae_metric(d, GT);
            Mmse(i,j) = mse_metric(d, GT);
        else
            Mmae(i,j) = NaN;
        end
    end
end

end