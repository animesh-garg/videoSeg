% Script to generate priors for our video segementatino algorithm
%   Uses the moseg_dataset and 3rd party NCuts code
%   Author: Jeff 
function [initialSegment, initialFlows, initialWeights, avgColorF, avgColorB] = ...
    generate_priors(videoName)

%% Get all images in a video sequence
startup;

videoPath = sprintf('%s%s', datasetPath, videoName);
folderContents = dir(videoPath);
videoFrames = {};
  
% read in all of the frames
k = 1;
for j = 1:size(folderContents,1)
    name = folderContents(j).name;
    if length(name) > 4 && ...
        strcmp(name((size(name,2)-3):size(name,2)), '.jpg')
        
        name = sprintf('%s/%s', videoPath, name);
        fprintf('Loading image %s ...\n',name);
        videoFrames{k} = imread(name);
        k = k+1;
    elseif length(name) > 6 && ...
        strcmp(name((size(name,2)-5):size(name,2)), '01.pgm')
        name = sprintf('%s/%s',videoPath, name);
        initialSegment = imread(name);
        
        % pick out only one of the possible object to segment, and store
        % initial segmentation as 0 or 1 valued array
        initialSegment(initialSegment < 255) = 0;      
        initialSegment(initialSegment==255) = 1;
    end
end
   
% downsample the images
[m, n, c] = size(videoFrames{1});
nFrames = size(videoFrames,2);
mDown = m / d;
nDown = n / d;
for i = 1:nFrames
    videoFrames{i} = imresize(videoFrames{i}, [mDown nDown]);
end
initialSegment = imresize(initialSegment, [mDown nDown]);

% get the initial flows
initialFlows = cell(1, nFrames-1);
optical_flow = vision.OpticalFlow('OutputValue', ...
        'Horizontal and vertical components in complex form', ...
        'ReferenceFrameSource', 'Input port', ...
        'Method', 'Lucas-Kanade');
for j =2:nFrames
    initialFlows{j-1} = step(optical_flow, ...
        double(rgb2gray(videoFrames{j})), ...
        double(rgb2gray(videoFrames{j-1})));
end

% convert flows to u,v and then to max likely weights
initialWeights = cell(1,size(initialFlows,2));
for j = 1:(nFrames-1)
    fprintf('Computing initial weights for frame %d...\n', j);
        
    U = real(initialFlows{j});
    V = imag(initialFlows{j});    
    initialWeights{j} = uv_to_weights(U, V, b);
end

%% Compute the appearance model constants
avgColorF = zeros(1,3);
avgColorB = zeros(1,3);

nFrames = size(videoFrames,2);
nPixFrame = mDown * nDown;
nPixForeground = sum(sum(initialSegment));
nPixBackground = nPixFrame - nPixForeground;
nPixTotal = nFrames * nPixFrame;

% get averages for appearance model cost
avgColorF(1) = (1.0 / nPixForeground) * ...
    sum(sum(videoFrames{1}(:,:,1).*initialSegment));
avgColorF(2) = (1.0 / nPixForeground) * ...
    sum(sum(videoFrames{1}(:,:,2).*initialSegment));
avgColorF(3) = (1.0 / nPixForeground) * ...
    sum(sum(videoFrames{1}(:,:,3).*initialSegment));

avgColorB(1) = (1.0 / nPixBackground) * ...
    sum(sum(videoFrames{1}(:,:,1).*(1-initialSegment)));
avgColorB(2) = (1.0 / nPixBackground) * ...
    sum(sum(videoFrames{1}(:,:,2).*(1-initialSegment)));
avgColorB(3) = (1.0 / nPixBackground) * ...
    sum(sum(videoFrames{1}(:,:,3).*(1-initialSegment)));

end
