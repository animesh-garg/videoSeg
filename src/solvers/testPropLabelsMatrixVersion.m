
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
startFrame = 15;
sequenceName = 'bird_of_paradise';
videoStruct = read_data(sequenceName, T, startFrame);
d = 1.0 / 20.0;

for t = 1:length(videoStruct.I)
    videoStruct.I{t} = imresize(videoStruct.I{t}, d);
end
videoStruct.X1 = imresize(videoStruct.X1, d);
videoStruct.X1 = videoStruct.X1(:,:,1);
%[videoStruct] = generateSyntheticDataMovement(Im,2,2,3,0.01);

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

% THINGS THAT WORK - tested on parachute with 3 frames
% [1 1e-2 1] tracks incorrectly
% [1e10 1e-2 1] tracks incorrectly

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
% set penalties
lambda = ones(1, 7);
lambda(1) = 1e0;
lambda(2) = 1e-1;

lambda(3) = 1e0;
lambda(4) = 1e2;
lambda(5) = 1e-1;
lambda(6) = 1e0;
lambda(7) = 0;

sigma = 0.8;

params.lambda =  lambda;
params.sigma = sigma;
params.window = window;
params.spatial_nbd_size = spatial_nbd_size;
iters = 1;

% loading params
useL2Penalty = false;

loadCSV = true;
saveCSV = false;

debug = true;


for k = 1:iters
    tic;
    X = solveLabels(X, U, V, videoStruct, T, lambda, spatial_nbd_size, ...
        window, useL2Penalty, loadCSV, saveCSV, debug);
    elapsed = toc;
    fprintf('Label solver took %f sec for iteration %d\n', elapsed, k);    

    tic;
    W = solveWeightsAvgMomentum(X, W, videoStruct, T, lambda, window, ...
        useL2Penalty, loadCSV, saveCSV, debug);
    elapsed = toc;
    fprintf('Weight solver took %f sec for iteration %d\n', elapsed, k);
    
    for t = 1:(T-1)
        [U{t},V{t}] = weights_to_uv(W{t});
    end
    
    loadCSV = true;
    saveCSV = false;
end


%% Get results

figure
I = videoStruct.I;

for t=1:T
    imgtmp = I{t};
    edgeTemp = edge(X{t});
    boundaryIdx = find (edgeTemp>0.5);
    [boundaryI, boundaryJ] = ind2sub(size(edgeTemp), boundaryIdx);
    C = size(I{t},3);    
    for p = 1:size(boundaryI)        
        imgtmp(boundaryI(p), boundaryJ(p),:) = 255;
    end
    subplot(2,T, t)
    imagesc(uint8(imgtmp));%hold on;pause(0.5);
    
    if (t<T)
        subplot (2,T,T+t)
        quiver(U{t},V{t});
        set (gca, 'YDIR', 'reverse');
    end
end


%%
figure;
subplot(1,2,1);
quiver(U{1},V{1});
set(gca,'YDir','reverse');
subplot(1,2,2);
imshow(videoStruct.I{1});

 