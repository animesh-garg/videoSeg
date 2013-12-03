function [U, V] = weights_to_uv(weights)
%WEIGHTS_TO_UV Summary of this function goes here
%   Detailed explanation goes here

    [m, n, win, win2] = size(weights);
    
    U = zeros(m, n);
    V = zeros(m, n);
    d = floor(win / 2);
    
    for i = 1:m
        for j = 1:n
            for a = 1:win
                for b = 1:win
                   U(i,j) = U(i,j) + weights(i,j,a,b) * (a-d-1);
                   V(i,j) = V(i,j) + weights(i,j,a,b) * (b-d-1);
                end
            end
        end
    end
    
    
end

