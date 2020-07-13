function d_new = quadratic_process(d, DSI, offset)
% Function which generates better estimate of disparity
% Uses quadratic interpolation, applied on DSI
%
% d is the disparity image
% DSI is the disparity space image
% 
% d_new is the interpolated disparity values

[A, B, C] = quad_int(DSI);
d_new = NaN(size(d));

for i = 1:size(A,1)
    for j = 1:size(A,2)
        a = A(i,j);
        b = B(i,j);
        c = C(i,j);
        
        % optimum of quadratic polynomial
        % also add offset
        if ~isnan(a)
            d_new(i,j) = -b / (2*a) + offset;
        else
            d_new(i,j) = d(i,j);
        end
    end
end