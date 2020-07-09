function [d, sim, DSI] = classic_method(L, R, d_range, H)
% Function for solving stereo matching using classic algorithm
% [d, sim, DSI] = classic_method(L, R, d_range, H)
%
% L, R are the image pair
% d_range = [d_min, d_max] is the range of feasible disparities
% H is the half-length of the correlation window

% ensure images are monochromatic
L = imono(L);
R = imono(R);

% get some parameters
dmin = d_range(1);
dmax = d_range(2);
n = size(L, 1);
m = size(L, 2);
D = dmax - dmin + 1;

% define return values
DSI = NaN(n, m, D);

for i = 1:n
    if mod(i,10) == 0
        fprintf('i = %d\n', i);
    end
    
    if i-H > 0 && i+H <= n % check validity of row index
        for j = 1:m
            if j-H > 0 && j+H <= m % validity of col index
                for dis = dmin:dmax
                    if j-dis-H > 0 % feasible matching
                        % match patches in L and R using ZNCC
                        Lp = L(i-H:i+H, j-H:j+H);
                        Rp = R(i-H:i+H, j-dis-H:j-dis+H);
                        DSI(i,j,dis-dmin+1) = my_zncc(Lp, Rp);
                    end
                end
            end
        end
    end
end

% get disparity estimates from disparity state image
[sim, d] = max(DSI, [], 3);

% if all values are NaN, set index of max to NaN
for i = 1:n
    for j=1:m
        if isnan(sim(i,j))
            d(i,j) = NaN;
        end
    end
end

% adjust for offset
d = d + (dmin - 1);

end