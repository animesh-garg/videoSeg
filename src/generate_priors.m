function [initU, initV, initW, avgColorF, avgColorB] = ...
    generate_priors(videoStruct)
% Takes the videoStruct as Input and 
% Outputs : initU, initV, initW, avgColorF, avgColorB
% initU, intiV : (T-1) X height X width flow field matrices
% initW : (T-1) X height X width X (2b+1) X (2b+1) weight vectors
% avgColorF, avgColorB : avg pixel color for foreground and background
% respectively

%% initializing some variable values
startup;
I = videoStruct.I;
X1 = videoStruct.X1;
%X1 = double(videoStruct.X1);
T = length(videoStruct.I);
[ydim,xdim,~] = size(videoStruct.I{1});
initU = cell(T-1,1);
initV =  cell(T-1,1);
initW = cell(T-1,1);

%% get the initial flows

for t = 1:T-1
    %compute flow using Horn-Schunk method and Flow function in Piotr's
    %toolbox
    fprintf('Computing flows for frame %d...\n', t); 
    %[initU{t},initV{t}] = opticalFlow(rgb2gray(I{t}),rgb2gray(I{t+1}),{'type','SD'});
    %[initU{t},initV{t}] = opticalFlow(rgb2gray(I{t}),rgb2gray(I{t}),{'type','LK'});
    UV = estimate_flow_interface(I{t},I{t+1},'classic+nl-fast');
    initU{t} = UV(:,:,1);
    initV{t} = UV(:,:,2);
end

%% get the initial weights

for t = 1:(T-1)
    fprintf('Computing initial weights for frame %d...\n', t);   
    initW{t} = uv_to_weights(initU{t}, initV{t}, b);
end

%% Compute the appearance model constants

avgColorF = zeros(1,3);
avgColorB = zeros(1,3);

nPixFrame = ydim * xdim;
nPixForeground = sum(sum(X1));
nPixBackground = nPixFrame - nPixForeground;

% get averages for appearance model cost
% frame1 = double(I{1});
% avgColorF(1) = (1.0 / nPixForeground) * ...
%     sum(sum(squeeze(frame1(:,:,1)).*X1));
% avgColorF(2) = (1.0 / nPixForeground) * ...
%     sum(sum(squeeze(frame1(:,:,2)).*X1));
% avgColorF(3) = (1.0 / nPixForeground) * ...
%     sum(sum(squeeze(frame1(:,:,3)).*X1));
% 
% avgColorB(1) = (1.0 / nPixBackground) * ...
%     sum(sum(squeeze(frame1(:,:,1)).*(1-X1)));
% avgColorB(2) = (1.0 / nPixBackground) * ...
%     sum(sum(squeeze(frame1(:,:,2)).*(1-X1)));
% avgColorB(3) = (1.0 / nPixBackground) * ...
%     sum(sum(squeeze(frame1(:,:,3)).*(1-X1)));

end
