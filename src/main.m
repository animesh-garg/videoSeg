% MAIN script for EE227B project

% set vars
iters = 10;
T = 5;
window = 2;
sequenceName = 'cars1';

% set penalties
lambda = ones(1, 7);
lambda(3) = 1e0;
lambda(4) = 1e0;
lambda(5) = 2e1;
lambda(6) = 1e1;
lambda(7) = 1e-2;
useL2Penalty = false;
debug = false;

video = read_data(sequenceName);
[initU, initV, initW, avgColorF, avgColorB] = generate_priors(video);

U = initU;
V = initV;

for k = 1:iters
    % some way to run label solver
    [U, V] = solveWeightsHornSchunk(X, U, V, video, T, lambda, window, ...
        useL2Penalty, debug);
end
