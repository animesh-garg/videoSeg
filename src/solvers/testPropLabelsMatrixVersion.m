
clear all, 
close all
clc

%% add paths
%addpath( './utils/')
startup();

%% create test data
%[videoStruct, T] = genSimpleTestVideo();
Im = imread('chess.jpg');
T = 3;
sequenceName = 'parachute';
videoStruct = read_data(sequenceName, T);

for t = 1:length(videoStruct.I)
    videoStruct.I{t} = imresize(videoStruct.I{t}, 1.0 / 16.0);
end
videoStruct.X1 = imresize(videoStruct.X1, 1.0 / 16.0);
videoStruct.X1 = videoStruct.X1(:,:,1);
%[videoStruct] = generateSyntheticDataMovement(Im,2,2,3,0.01);

%T = 5

%% get flows
% set vars
window = 2; %value of b temp nbd
spatial_nbd_size = 1; %spatial nbd
[m, n, c] = size(videoStruct.I{1});
T = length(videoStruct.I);
T = 3;

% THINGS THAT WORK - tested on square with 3 frames
% [0 1e-2 1] tracks square
% [1 1e-2 1] tracks square

% set penalties
lambda = ones(1, 7);
lambda(1) = 1e0;
lambda(2) = 1e-2;
lambda(3) = 1e0;

lambda(4) = 1e0;
lambda(5) = 1e3;
lambda(6) = 1e1;
lambda(7) = 0;

sigma = 0.8;
useL2Penalty = false;
loadCSV = false;
saveCSV = true;
debug = true;

[initU, initV, initW, avgColorF, avgColorB] = generate_priors(videoStruct);

% video = read_data(sequenceName);
% [initU, initV, initW, avgColorF, avgColorB] = generate_priors(video);

U = initU;
V = initV;
W = initW;
initX = cell(T, 1);
initX{1} = videoStruct.X1;
for t = 2:T
    initX{t} = zeros(m,n);
end
X = initX;

%% define params
params.lambda =  lambda;
params.sigma = sigma;
params.window = window;
params.spatial_nbd_size = spatial_nbd_size;
iters = 1;
 
for k = 1:iters
    tic;
    X = solveLabels(X, U, V, videoStruct, T, lambda, spatial_nbd_size, ...
        window, useL2Penalty, loadCSV, saveCSV, debug);
    elapsed = toc;
    fprintf('Label solver took %f sec for iteration %d\n', elapsed, k);
    
    tic;
    [U, V] = solveWeightsHornSchunk(X, U, V, videoStruct, T, lambda, window, ...
        useL2Penalty, loadCSV, saveCSV, debug);
    elapsed = toc;
    fprintf('Flow solver took %f sec for iteration %d\n', elapsed, k);
    
%     loadCSV = true;
%     saveCSV = false;
end
%%
visualizeSegmentationResult(videoStruct, X);
 