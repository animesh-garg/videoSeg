function [] = visualizeSegmentationResult(vidStruct,X)
%VISUALIZESEGMENTATIONRESULT Visualize obtained segmentations.
%   Plays the video by dimming the pixels out of the labelled foreground

T = length(X);
I = vidStruct.I;
brightnessFraction = 0.2;

for t=1:T
    imgtmp = zeros(size(I{t}));
    C = size(I{t},3);
    for c=1:C
        imgtmp(:,:,c) = uint8(I{t}(:,:,c).*uint8(X{t}>0.5) + brightnessFractions*(I{t}(:,:,c).*uint8(X{t}<=0.5)));
    end
    imagesc(uint8(imgtmp));hold on;pause(0.2);
end

end

