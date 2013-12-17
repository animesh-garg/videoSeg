% junk script to test flow computation

% setup the problem
%video = genSimpleTestVideo();
clear;
clc;
close all;
times = [];

% recording times for experiment w/ 16x16 image of downsampled chessboard
% versus time win=1, 16x16
% versus weight window size 16x16, T=4
% versus image size win=1, T=4

%%
T = 5;
image = imread('utils/chess.jpg');
video = generateSyntheticDataMovement(image, -1*sin(1:T), 1*cos(1:T), T, 0.01);

figure;
for i = 1:T
    subplot(1,T,i);
    imshow(video.I{i});
end

[m, n, c] = size(video.I{1});


%% calculate new weights
U_init = cell(1,T);
V_init = cell(1,T);
W_init = cell(1,T);

windowSize = 1;

for t = 1:T
    U_init{t} = zeros(m, n);
    V_init{t} = zeros(m, n);
    W_init{t} = uv_to_weights(U_init{t}, V_init{t}, windowSize);
end

%% SOLVE

% set penalties
lambda = ones(1, 7);
lambda(3) = 1e0;
lambda(4) = 1e0;
lambda(5) = 2e1;
lambda(6) = 1e1;
lambda(7) = 1e-2;

useL2Penalty = false;
loadCSV = false;
saveCSV = true;
debug = false;

% solve
tic;
[newW] = solveWeightsHornSchunk(video.gTruth.X, U_init, V_init, ...
    video, T, lambda, windowSize, useL2Penalty, loadCSV, saveCSV, debug);
elapsed = toc;

times = [times elapsed];

fprintf('Optimization took %f sec\n', elapsed); 

%% visualize
figure;
for i = 1:(T-1)
    %figure;
    subplot(1,(T-1),i);
    [newU{i}, newV{i}] = weights_to_uv(newW{i});
    visualizeFlow(video.I{i}, newU{i}, newV{i});
end

%% plot
figure;
plot([4 8 16 24], times);
title('Time to Solve Versus Image Size', 'FontSize', 20);
xlabel('Image Width (pixels)', 'FontSize', 15);
ylabel('Time to Solve (sec)', 'FontSize', 15);