function weights = uv_to_weights(U, V, b)
%UV_TO_WEIGHTS Summary of this function goes here
%   Detailed explanation goes here
%TODO - use interpolation instead of rounding
% b: half size of the window
    [m, n] = size(U);

    U(U >  b) =  b;
    U(U < -b) = -b;
    V(V >  b) =  b;
    V(V < -b) = -b;
    U = int8(round(U));
    V = int8(round(V));
        
    weights = zeros(m, n, 2*b+1, 2*b+1);
    for x = 1:m
        for y = 1:n
            weights(x,y,U(x,y)+b+1,V(x,y)+b+1) = 1;
        end
    end
end

