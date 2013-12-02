function [cost] = avgValueBasedCost(vid)
%AVGVALUEBASEDCOST Gives foreground cost for each pixel
%   Input : VideoStruct
%   Output : cell of length T. Each cell element is a height X widhth
%   matrix with value = foreground cost for that pixel (between [0 1]). The
%   background cost is 1 - foreground Cost
%   Even though this is very fast, I would suggest use of the randomForest based model over this 
%   as it is more reliable as it handles variations within foreground and background well

%% Initializing Variables
I = vid.I;
X1 = double(vid.X1);
T = length(I);
[ydim,xdim,~] = size(I{1});

%% Compute the appearance model constants

avgColorF = zeros(1,3);
avgColorB = zeros(1,3);

nPixFrame = ydim * xdim;
nPixForeground = sum(sum(X1));
nPixBackground = nPixFrame - nPixForeground;

% get averages for appearance model cost
frame1 = double(I{1});
avgColorF(1) = (1.0 / nPixForeground) * ...
    sum(sum(squeeze(frame1(:,:,1)).*X1));
avgColorF(2) = (1.0 / nPixForeground) * ...
    sum(sum(squeeze(frame1(:,:,2)).*X1));
avgColorF(3) = (1.0 / nPixForeground) * ...
    sum(sum(squeeze(frame1(:,:,3)).*X1));

avgColorB(1) = (1.0 / nPixBackground) * ...
    sum(sum(squeeze(frame1(:,:,1)).*(1-X1)));
avgColorB(2) = (1.0 / nPixBackground) * ...
    sum(sum(squeeze(frame1(:,:,2)).*(1-X1)));
avgColorB(3) = (1.0 / nPixBackground) * ...
    sum(sum(squeeze(frame1(:,:,3)).*(1-X1)));

%% Compute Costs
cost = cell(T,1);
cost{1} = ones(ydim,xdim)-X1;
for t=2:T
    frame = double(I{t});
    tmpF(1,1,:) = avgColorF;
    diffFg = frame - repmat(tmpF,[ydim xdim 1]);
    distFg = sum(abs(diffFg),3);
    
    tmpB(1,1,:) = avgColorB;
    diffBg = frame - repmat(tmpB,[ydim xdim 1]);
    distBg = sum(abs(diffBg),3);
    
    cost{t} = distFg./(distFg+distBg);
end

end

