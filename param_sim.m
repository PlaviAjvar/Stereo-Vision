function [Mlam, Msat, Mmae] = param_sim(image, method)
% Function used to simulate parametrized algorithms
% [Mlam, Msat, Mmae] = param_sim(image, method)
%
% method = {'ScanlineDP', 'SGM', 'LoopyBP', 'MultiscaleBP'}
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
        if Mlam(i,j) < Msat(i,j)
            DSI = fast_stereo(L, R, H, [dmin, dmax], method, [Mlam(i,j), Msat(i,j)]);
            [sim, d] = min(DSI, [], 3);
            d(isnan(sim)) = NaN;
            d = d + offset;
            d = d + (dmin - 1);
            Mmae(i,j) = dif_metric(d, GT);
        else
            Mmae(i,j) = NaN;
        end
    end
end

end