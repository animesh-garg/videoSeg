function [Cost] = getUnaryPotential(imSequence, initLabel)
%Input: 
%   init_label: the matrix containing the label for initial image (initial 
%       segmentation) 
%   imSequence: sequence of image of length T
%   T: number of time steps in future for prediction. (parameter)


[M, N, T] = size (imSequence);

%init_label seeds the similarity for pixel intensities defined to be
%foregraound

%find indices for foreground labels
[Idx, ~] = find(initLabel);

%get foreground pixels in initial image in the sequence. 
initImage = reshape(imSequence(:,:, 1), M,N);
pixelsInFG = initImage(Idx);
sizeFG = length(pixelsInFG);

% For every pixel intensity calculate similarity
for t = 1: T
   for i= 1:M
    for j = 1:N
        %assuming Grayscale for now
        pixelDist = imSequence (i,j,t)*ones(sizeFG,1) - pixelsInFG ; 
        if (min(pixelDist) < threshEPS))
            %if similar to FG assign cost as min distance (normalized)
            Cost (i,j,t) = min(pixelDist)/sum(pixelDist); 
        else
            %if similar to FG assign a large cost
            Cost (i,j,t) = 1;
        end
    end
   end   
end


% Similarity can be pixel wise or patch wise.
% facier ways would be withint some mean of a guassian mixture

end

