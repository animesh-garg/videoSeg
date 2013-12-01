% junk script to test flow computation

% setup the problem
video = genSimpleTestVideo();
figure;
subplot(1,3,1);
imshow(video.I{1});
subplot(1,3,2);
imshow(video.I{2});
subplot(1,3,3);
imshow(video.I{3});

d = 1;
[m, n, c] = size(video.I{1});
T = video.T;
X = cell(1,T);
X{1} = video.X1;
X{2} = zeros(m, n);
X{2}(3:4,3:4) = ones(2,2);
X{3} = zeros(m, n);
X{3}(4:5,4:5) = ones(2,2);

%% calculate new weights
W = cell(1,T);
W{1} = zeros(m, n, 2*d+1, 2*d+1);
W{2} = zeros(m, n, 2*d+1, 2*d+1);
for i = 1:m
    for j = 1:n
        W{1}(i,j,2,2) = 1;
    end
end

lambda = ones(1, 6);
newW = solveWeights(X, W, video, lambda, d);

%% visualize
figure;
[U, V] = weights_to_uv(newW{1});
visualizeFlow(video.I{1}, U, V);

figure;
[U, V] = weights_to_uv(newW{2});
visualizeFlow(video.I{2}, U, V);