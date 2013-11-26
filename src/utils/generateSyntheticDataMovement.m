function [videoStruct] = generateSyntheticDataMovement(seedImg,yVel,xVel,T,resizeFactor)
%GENERATESYNTHETICDATA Gives a videoStruct augmented with ground truth segmentations
%   INPUTS -
%   seedImg : A simple image which corresponds to the moving foreground
%   yVel, xVel : velocity with which you want the image to move. Can be
%   array of velocities for each time T or a sigle entry for a constant
%   velocity
%   T : number of desired frames
%   resizeFactor : optional argument for scaling down image
%
% OUTPUTS - 
%   creates a video for an image 3 times the size of the seedImage*resizeFactor and
%   returns the corresponding videoStruct
%   example : generateSyntheticDataMovement(seedImg,2+2*sin(1:20),2+2*cos(1:20),20,0.2);


%% Initializing the frame 1 image
%seedImg = rgb2gray(seedImg);
[ydim,xdim,C] = size(seedImg);
initFrame = zeros(3*ydim,3*xdim,C);
X1 = zeros(3*ydim,3*xdim);
initFrame(ydim+1:2*ydim,xdim+1:2*xdim,:) = seedImg;
X1(ydim+1:2*ydim,xdim+1:2*xdim)=1;
videoStruct = struct;

%% Generating velocity matrix


xVel = round(xVel);
yVel = round(yVel);

if(size(yVel,1)>1)
    yVel = yVel';
    xVel = xVel';
end

if(size(yVel,2) ~= size(xVel,2))
    disp('Error! xVel and yVel not of same sizes');
    return;
end

if(size(yVel,2)==1)
    yVel = repmat(yVel,1,T);
    xVel = repmat(xVel,1,T);
elseif(size(yVel,2)<(T-1))
    disp('Error! xVel and yVel not enough for T frames');
    return;
end

%% Resizing
if(nargin < 5)
    resizeFactor = 0.3;
end

I_t = double(initFrame);
I_t = imresize(I_t,resizeFactor);
X1 = round(imresize(X1,resizeFactor));

%% Loop to generate video synthetic data

I = cell(T,1);
for t=1:T
    
    w = max(abs([yVel(t) xVel(t)]));
    convFilter = zeros(2*w+1,2*w+1);
    convFilter(yVel(t)+w+1,xVel(t)+w+1)=1;
    
    imagesc(uint8(I_t));
    pause(0.1);
    
    I{t} = uint8(I_t);
    for c = 1:C
        I_t(:,:,c) = conv2(I_t(:,:,c),convFilter,'same');
    end
end

videoStruct.I = I;
videoStruct.X1 = X1;

%% Computing ground truth segmentations and flows for foreground
U = cell(T-1,1);
V = cell(T-1,1);
X = cell(T,1);

X_t = X1;
X{1}=X1;

for t=1:T-1
    
    w = max(abs([yVel(t) xVel(t)]));
    convFilter = zeros(2*w+1,2*w+1);
    convFilter(yVel(t)+w+1,xVel(t)+w+1)=1;
    
    U{t} = zeros(size(X_t));
    V{t} = zeros(size(X_t));
    U{t}(X_t > 0) = xVel(t);
    V{t}(X_t > 0) = yVel(t);
    
    X_t = conv2(X_t,convFilter,'same');
    X{t+1} = X_t;
end

gTruth.X = X;
gTruth.U = U;
gTruth.V = V;

videoStruct.gTruth = gTruth;

end