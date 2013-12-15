
clear all, 
close all
clc

%% add paths
%addpath( './utils/')
startup();

%% create test data
%[videoStruct, T] = genSimpleTestVideo();
Im = imread('chess.jpg');
[videoStruct] = generateSyntheticDataMovement(Im,1,1,5,0.001);

%T = 5;

%% get flows
% set vars
iters = 10;
window = 2; %value of b temp nbd
spatial_nbd_size = 1; %spatial nbd
sequenceName = 'cars1';
T = length(videoStruct.I);

% set penalties
lambda = ones(1, 7);
lambda(3) = 1e0;
lambda(4) = 1e0;
lambda(5) = 2e1;
lambda(6) = 1e1;
lambda(7) = 1e-2;
sigma = 0.8;
useL2Penalty = false;
debug = false;

[initU, initV, initW, avgColorF, avgColorB] = generate_priors(videoStruct);

% video = read_data(sequenceName);
% [initU, initV, initW, avgColorF, avgColorB] = generate_priors(video);

U = initU;
V = initV;
W = initW;

%% define params
params.lambda =  lambda;
params.sigma = sigma;
params.window = window;
params.spatial_nbd_size = spatial_nbd_size;

% [U, V] = solveWeightsHornSchunk(X, U, V, video, T, lambda, window, ...
%         useL2Penalty, debug);

%[U, V] = solveWeightsHornSchunk(X, U, V, videoStruct, T, lambda, window, ...
%       useL2Penalty, debug);

    
%W = uv_to_weights(U,V,window);
 
X = propagate_labels(videoStruct, W, params);

visualizeSegmentationResult(videoStruct, X);
 
%% propagate labels