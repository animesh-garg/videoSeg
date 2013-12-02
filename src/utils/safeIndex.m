function di = safeIndex(i, d, M)
%SAFEINDEX Compute safe index of i+d, snapping to image bounds
    di = min(max(i+d,1), M);
end

