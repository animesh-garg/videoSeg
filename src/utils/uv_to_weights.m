function weights = uv_to_weights(U, V, win)
%UV_TO_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here
%TODO - use interpolation instead of rounding
% b: half size of the window
    [m, n] = size(U);

    U(U >  win) =  win;
    U(U < -win) = -win;
    V(V >  win) =  win;
    V(V < -win) = -win;
        
    weights = zeros(m, n, 2*win+1, 2*win+1);
    for x = 1:m
        for y = 1:n
            a_upper = ceil(U(x,y));
            a_lower = floor(U(x,y));
            b_upper = ceil(V(x,y));
            b_lower = floor(V(x,y));
            
            u_upper = (U(x,y) - a_lower) / (a_upper - a_lower);
            v_upper = (V(x,y) - b_lower) / (b_upper - b_lower);
            u_lower = (a_upper - U(x,y)) / (a_upper - a_lower);
            v_lower = (b_upper - V(x,y)) / (b_upper - b_lower);
            
            if a_upper == a_lower
                u_lower = 1;
                u_upper = 1;
            end    
            if b_upper == b_lower 
                v_lower = 1;
                v_upper = 1;
            end
            weights(x,y,a_lower+win+1,b_lower+win+1) = ...
                u_lower * v_lower;
            weights(x,y,a_upper+win+1,b_lower+win+1) = ...
                u_upper * v_lower;
            weights(x,y,a_lower+win+1,b_upper+win+1) = ...
                u_lower * v_upper;
            weights(x,y,a_upper+win+1,b_upper+win+1) = ...
                u_upper * v_upper;
        end
    end
end

