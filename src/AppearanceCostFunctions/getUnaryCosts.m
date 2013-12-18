function [Cost] = getUnaryPotential(imSequence, initLabel, threshEPS, pixelWise)
%Input: 
%   init_label: the matrix containing the label for initial image (initial 
%       segmentation) 
%   imSequence: sequence of image of length T
%   T: number of time steps in future for prediction. (parameter)
%   pixelWise: Optional argument (0 or 1)
%              Decides whether to compute the similarities pixelwise (1) or 
%              patchwise (0). Default patchSize 5x5



[M, N, T] = size (imSequence);

if ~exist('pixelWise')
    pixelWise = 1;
end

patchSize = 2; %this creates a 5 x 5 patch
nbdSize = (2*patchSize+1)^2;
%init_label seeds the similarity for pixel intensities defined to be
%foregraound

%find indices for foreground labels

if (pixelWise ==1)
    [Idx, ~] = find(initLabel);
    sizeFG = length(Idx);
    fprintf('Calculating Unary Potentials Pixel Wise \n')
    %get foreground pixels in initial image in the sequence. 
    initImage = reshape(imSequence(:,:, 1), M,N);
    pixelsInFG = initImage(Idx);
    
else
    fprintf('Calculating Unary Potentials Patch Wise \n')
    [I,J] = find(initLabel);
    sizeFG = length(I);
    initImage = reshape(imSequence(:,:, 1), M,N);
    
    for p = 1:sizeFG
        %imageIDX = ind2sub([M,N], idx(p));
        for i = (-patchSize): patchSize
            for j = (-patchSize): patchSize
                nbd_i= min(max(i+I(p),1), M);
                nbd_j = min(max(j+J(p),1), N);
                %GRAYSCALE ASSUMPTION - TODO: Color
                pixelsInFG (p) = pixelsInFG(p) + ...
                    initImage(nbd_i, nbd_j)/nbdSize;                
            end
        end
    end
end

% For every pixel intensity calculate similarity
for t = 1: T
   for i= 1:M
    for j = 1:N
        %GRAYSCALE ASSUMPTION - TODO: Color
        pixelDist = imSequence (i,j,t)*ones(sizeFG,1) - pixelsInFG ; 
        if (min(pixelDist) < threshEPS)
            %if similar to FG assign cost as min distance (normalized)
            Cost (t,i,j) = min(pixelDist)/sum(pixelDist); 
        else
            %if similar to FG assign a large cost
            Cost (t,i,j) = 1;
        end
    end
   end   
end


%TODO % Similarity can be pixel wise or patch wise. 
% facier ways would be withint some mean of a guassian mixture

end

