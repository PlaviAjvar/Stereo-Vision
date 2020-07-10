function [Mlam, Msat, Mmae] = DP_param_sim(image)
% Function used to simulate scanline DP with smoothness constraint
% [Mlam, Msat, Mmae] = DP_param_sim(image)
%
% Mlam, Msat are the matrices of scalar/saturation pairs (meshgrid)
% Mmae is the mean absolute error for GT and d obtained, using given pair

lam  = [0.05 0.2 0.4 1];
vsat = [0.1  0.5 1   2];

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
fileID = fopen(strcat(image,'\dmin.txt'));
offset = fscanf(fileID, '%d') / 2;
GT = GT + offset;
    
for i = 1:size(lam,2)
    for j = 1:size(vsat,2)
        DSI = fast_stereo(L, R, H, [dmin, dmax], 'SmoothDP', [Mlam(i,j), Msat(i,j)]);
        [sim, d] = min(DSI, [], 3);
        d(isnan(sim)) = NaN;
        d = d + offset;
        d = d + (dmin - 1);
        Mmae(i,j) = dif_metric(d, GT);
    end
end

end