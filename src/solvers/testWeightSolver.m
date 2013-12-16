% junk script to test flow computation

% setup the problem
%video = genSimpleTestVideo();

T = 5;
image = imread('src/utils/chess.jpg');
video = generateSyntheticDataMovement(image, -2*sin(1:T), 2*cos(1:T), T, 0.01);

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

windowSize = 2;

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
loadCSV = true;
saveCSV = false;
debug = false;

% solve
tic;
[newW] = solveWeightsAvgMomentum(video.gTruth.X, W_init, ...
    video, T, lambda, windowSize, useL2Penalty, loadCSV, saveCSV, debug);
elapsed = toc;

fprintf('Optimization took %f sec\n', elapsed); 

%% visualize
figure;
for i = 1:(T-1)
    %figure;
    subplot(2,(T-1)/2,i);
    [newU{i}, newV{i}] = weights_to_uv(newW{i});
    visualizeFlow(video.I{i}, newU{i}, newV{i});
end