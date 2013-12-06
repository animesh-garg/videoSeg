% junk script to test flow computation

% setup the problem
%video = genSimpleTestVideo();

T = 3;
image = imread('src/utils/chess.jpg');
video = generateSyntheticDataMovement(image, [2, 0, 0], [2, -2, 0], T, 0.05);

figure;
for i = 1:T
subplot(1,T,i);
imshow(video.I{i});
end

d = 2;
[m, n, c] = size(video.I{1});

%% calculate new weights
U_init = cell(1,T);
V_init = cell(1,T);
for t = 1:T
    U_init{t} = zeros(m, n);
    V_init{t} = zeros(m, n);
end

lambda = ones(1, 6);
tic;
newW = solveWeightsHornSchunk(video.gTruth.X, U_init, V_init, video, T, lambda, d);
elapsed = toc;

fprintf('Optimization took %f sec\n', elapsed); 

%% visualize
figure;
[U, V] = weights_to_uv(newW{1});
visualizeFlow(video.I{1}, U, V);

figure;
[U, V] = weights_to_uv(newW{2});
visualizeFlow(video.I{2}, U, V);
%%
figure;
[U, V] = weights_to_uv(newW{3});
visualizeFlow(video.I{3}, U, V);