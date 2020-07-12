function [A, B, C] = quad_int(DSI)
% Function which interpolates quadratic through neighbor pixels
% [A, B, C] = quad_int(DSI)
%
% DSI is the disparity space image
%
% A, B, C are matrices which contain the pixel coefficients

A = NaN([size(DSI,1), size(DSI,2)]);
B = NaN([size(DSI,1), size(DSI,2)]);
C = NaN([size(DSI,1), size(DSI,2)]);

for y = 1:size(DSI,1)
    for x = 1:size(DSI,2)
        [simv, dv] = min(DSI(y,x,:));
        
        if dv ~= 1 && dv ~= size(DSI,3)
            sim_L = DSI(y,x,dv-1);
            sim_R = DSI(y,x,dv+1);
            
            % interpolation matrix form
            mat = [(dv-1)^2 (dv-1) 1; dv^2 dv 1; (dv+1)^2 dv+1 1];
            rhs = [sim_L; simv; sim_R];
            
            % solve for parameters
            coef = mat^(-1) * rhs;
            A(y,x) = coef(1);
            B(y,x) = coef(2);
            C(y,x) = coef(3);
        end
    end
end

end