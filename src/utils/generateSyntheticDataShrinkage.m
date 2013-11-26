function [videoStruct] = generateSyntheticDataShrinkage(seedImg, T, resizeFactor, shrinkFactor)
%GENERATESYNTHETICDATASHRINKAGE Gives a videoStruct augmented with ground truth segmentations
%   INPUTS -
%   seedImg : A simple image which corresponds to the moving foreground
%   T : number of desired frames
%   resizeFactor : optional argument for scaling down image
%   shrinkFactor : fraction by which the seed image shrinks in every
%   frame
%
% OUTPUTS - 
%   creates a video for an image the size of the seedImage*resizeFactor and
%   returns the corresponding videoStruct
%   example : generateSyntheticDataShrinkage(seedImg,20,0.8,0.95);


%% Initial frame
[ydim,xdim,C] = size(seedImg);
X1 = ones(ydim,xdim);
videoStruct = struct;
seedImg = imresize(seedImg,resizeFactor);
X1 = round(imresize(X1,resizeFactor));

%% Loop to generate image frames, segmentations and compute flows
I = cell(T,1);
U = cell(T,1);
V = cell(T,1);
X = cell(T,1);
I{1} = uint8(seedImg);
X{1} = X1;

for t=1:T-1
    %% Initializing to zeros
    imagesc(I{t});hold on;pause(0.2);
    I{t+1} = zeros(ydim,xdim,C);
    X{t+1} = zeros(ydim,xdim);
    U{t} = zeros(ydim,xdim);
    V{t} = zeros(ydim,xdim);
    
    [ydimOld,xdimOld,~] = size(seedImg);
    seedImg = imresize(seedImg,shrinkFactor);
    [ydimNew,xdimNew,~] = size(seedImg);
    
    %% Assigning segmentation mask for new frame
    X{t+1}(1:ydimNew,1:xdimNew)=ones(ydimNew,xdimNew);
    
    %% Pixel values for new frame
    I{t+1}(1:ydimNew,1:xdimNew,:) = seedImg;
    I{t+1} = uint8(I{t+1});
    
    %% Computing Flows
    deltaX = 1-xdimNew/xdimOld;
    deltaY = 1-ydimNew/ydimOld;
    U{t}(1:ydimOld,1:xdimOld) = -deltaX*repmat(1:xdimOld,ydimOld,1);
    V{t}(1:ydimOld,1:xdimOld) = -deltaY*repmat((1:ydimOld)',1,xdimOld);
    
end

%% Returning Video struct
gTruth.X = X;
gTruth.U = U;
gTruth.V = V;

videoStruct.gTruth = gTruth;
videoStruct.I = I;
videoStruct.X1 = X{1};
end

