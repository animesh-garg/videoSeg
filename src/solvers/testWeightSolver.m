% junk script to test flow computation

% setup the problem
%video = genSimpleTestVideo();

T = 3;
image = imread('src/utils/chess.jpg');
video = generateSyntheticDataMovement(image, [2, 0, 0], [2, -2, 0], T, 0.01);

figure;
for i = 1:T
subplot(1,T,i);
imshow(video.I{i});
end

d = 2;
% [m, n, c] = size(video.I{1});
% X = cell(1,T);
% X{1} = video.X1;
% X{2} = zeros(m, n);
% X{2}(3:m,3:n) = video.X1(1:(m-2),1:(n-2));
% X{3} = zeros(m, n);
% X{3}(5:m,5:n) = video.X1(1:(m-4),1:(n-4));

%% calculate new weights
W = cell(1,T);
for t = 1:T
    W{t} = zeros(m, n, 2*d+1, 2*d+1);
    for i = 1:m
        for j = 1:n
            W{t}(i,j,2,2) = 1;
        end
    end
end

lambda = ones(1, 6);
tic;
newW = solveWeightsNoMomentum(video.gTruth.X, W, video, T, lambda, d);
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