function [newW] = solveWeightsAvgMomentum(X, W, video, T, lambda, ...
    windowSize, useL2Penalty, loadCSV, saveCSV, debug)
%SOLVE_WEIGHT_AVG_MOMENTUM Solve for the flow weights given a pixel labeling (using
% cplex) with average momentum constraint
% The variable we end up optimizing is W_hat = [W Y Z R S]', where
% Y puts the U spatial constraints into epigraph form,
% Z puts the V spatial constraints into epigraph form,
% R puts the U momentum constraints into epigraph form, and
% Z puts the V momentum constraints into epigraph form
%   INPUTS -
%   X : tentative labels of pixels in the video sequence
%   W : initial flow weights to seed the optimization
%   video : struct containing the video data
%   lambda : vector of constant weights in the objective function
%   d : size of window
%
% OUTPUTS - 
%   returns the weights found by cplex
%   example : generateSyntheticDataMovement(seedImg,2+2*sin(1:20),2+2*cos(1:20),20,0.2);

    win = 2*windowSize+1; % full window dimension
    [M, N] = size(X{1});
    c = size(video.I{1}, 3);
    
    numWeightsPerPixel = win^2;
    numPixelsPerFrame = M * N;
    numPixels = (T-1) * numPixelsPerFrame; % no weights for last frame
    numWeightsPerFrame = numWeightsPerPixel * numPixelsPerFrame;
    numWeights = numPixels * numWeightsPerPixel;
    numMomentum = (T-2) * numPixelsPerFrame;
    numVariables = numWeights + numPixels + numPixels + numMomentum + numMomentum;
    numVariableTypes = 5;
    
    % Compute the average intensity for each of the images
    K = zeros(c,1);
    for i = 1:c
        for t = 1:T
           K(i) = ((t-1)/t) * K(i) + (1.0/t) * mean(mean(video.I{t}(:,:,i)));
        end
    end
    
    % Normalize lambdas
    lambda(3) = double(lambda(3)) / double((win^2)*M*N*(T-1));
    lambda(4) = double(lambda(4)) / double((win^2)*M*N*(T-1)*(sum(K)));
    lambda(5) = double(lambda(5)) / double(2*M*N*(T-1)); 
    lambda(6) = double(lambda(6)) / double(2*M*N*(T-2)); 
    
    % Compute all of the start indices
    startIndices = zeros(numVariableTypes,1);
    startIndices(1) = 1;                            % W
    startIndices(2) = startIndices(1) + numWeights; % Y
    startIndices(3) = startIndices(2) + numPixels;  % Z
    startIndices(4) = startIndices(3) + numPixels;   % R
    startIndices(5) = startIndices(4) + numMomentum; % S
    
    % T and F costs
    disp('Creating linear costs');
    P = zeros(numVariables,1);
    Q = zeros(numVariables,1);

    % create label difference P, color difference Q matrices
    index = 1;
    for t = 1:(T-1)
        X_t = X{t};
        X_tplus = X{t+1};
        I_t = double(video.I{t});
        I_tplus = double(video.I{t+1});
        for i = 1:M
            for j = 1:N
                for b = (-windowSize):windowSize
                    for a = (-windowSize):windowSize
                        i_tplus = safeIndex(i, a, M);
                        j_tplus = safeIndex(j, b, N);

                        P(index) = abs(X_t(i,j) - X_tplus(i_tplus,j_tplus));
                        for k = 1:c
                            Q(index) = Q(index) + ...
                                abs(I_t(i,j,1) - I_tplus(i_tplus,j_tplus,k));
                        end
                        index = index+1;
                    end
                end 
            end     
        end
    end
    
    P = sparse(P);
    Q = sparse(Q);
    
    if loadCSV
        C = load('weight_cache/C.mat');
        C = C.C;
        H = load('weight_cache/H.mat');
        H = H.H;
    else
        % C cost - create vector to sum auxiliary flow difference variables
        C = zeros(numVariables,1);
        C(startIndices(2):startIndices(4)-1) = ones(2*numPixels,1);
        C = sparse(C);

        % M cost - create vector to sum auxiliary momentum variables
        H = zeros(numVariables,1);
        H(startIndices(4):numVariables) = ones(2*numMomentum,1);
        H = sparse(H);
    end
    
    P = lambda(3)*P + lambda(4)*Q + lambda(5)*C + lambda(6)*H; % final linear cost term
    
    if saveCSV
        save('weight_cache/C.mat', 'C');
        save('weight_cache/H.mat', 'H');
    end
    clear Q;
    clear C;
    clear H;
    
    % Equality constraints - ensure weights sum to 1
    disp('Generating equality constraints');
    if loadCSV
        Aeq = load('weight_cache/Aeq.mat');
        Aeq = Aeq.Aeq;
        beq = load('weight_cache/beq.mat');
        beq = beq.beq;
    else
        Aeq = zeros(numPixels, numVariables);
        startIndex = 1;
        endIndex = startIndex + numWeightsPerPixel -1;
        for i = 1:numPixels
           Aeq(i,startIndex:endIndex) = ones(1, numWeightsPerPixel);
           startIndex = endIndex + 1;
           endIndex = startIndex + numWeightsPerPixel - 1;
        end
        beq = ones(numPixels, 1);
        Aeq = sparse(Aeq);
    end
    
    % Inequality constraints
    % 1. Weights greater than 0
    % 2. Y auxiliary variables are greater that +/- the U flow continuity
    % 3. Z auxiliary variables are greater that +/- the V flow continuity
    % 4. R auxiliary variables are greater that +/- the U avg flow diff
    % between frames
    % 5. S auxiliary variables are greater that +/- the V avg flow diff
    % between frames
    disp('Generating inequality constraints');
    
    if loadCSV
       Aineq = load('weight_cache/Aineq.mat');
       Aineq = Aineq.Aineq;
       bineq = load('weight_cache/bineq.mat');
       bineq = bineq.bineq;
    else
        numInequalityConstraints = numWeights + 2*numPixels + 2*numPixels + ...
        2*numMomentum + 2*numMomentum;
        Aineq = zeros(numInequalityConstraints, numVariables);
        bineq = zeros(numInequalityConstraints, 1);

        % Constraint 1
        disp('Generating inequality constraint 1');
        Aineq(1:numWeights,1:numWeights) = -eye(numWeights);
        index = numWeights+1;

        % helper vectors for computing U, V
        alpha = repmat(linspace(-windowSize, windowSize, win), [1, win]);
        beta = reshape(repmat(linspace(-windowSize, windowSize, win), [win, 1]), [1 win^2]);

        % Constraint 2 - U spatial consistency constraint
        disp('Generating inequality constraint 2');
        startIndex = 1;
        endIndex = startIndex + numWeightsPerPixel - 1;
        for t = 0:(T-2)
            for i = 1:M
                for j = 1:N
                    c = zeros(1, numVariables);
                    c(startIndex:endIndex) = alpha;

                    % generate the mean computing matrices
                    startX = safeIndex(j, -windowSize, N);
                    endX = safeIndex(j, windowSize, N);
                    startY = safeIndex(i, -windowSize, M);
                    endY = safeIndex(i, windowSize, M);
                    xInterval = endX-startX+1;
                    yInterval = endY-startY+1;
                    nSummed = xInterval*yInterval;
                    startX = numWeightsPerPixel*(startX-1)+1;
                    endX = numWeightsPerPixel*endX;

                    f = zeros(M, numWeightsPerPixel*N);
                    f(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                        repmat(alpha, [yInterval xInterval]);
                    f = reshape(f', [1 numWeightsPerFrame]);
                    g = zeros(1, numVariables);
                    g((t*numWeightsPerFrame+1):((t+1)*numWeightsPerFrame)) = f;

                    Y_index = numWeights + (j-1)+(i-1)*N+t*M*N + 1;
                    Aineq(index+0, :) = c - g;
                    Aineq(index+0, Y_index) = -1;

                    Aineq(index+1, :) = g - c;
                    Aineq(index+1, Y_index) = -1;

                    startIndex = endIndex + 1;
                    endIndex = startIndex + numWeightsPerPixel - 1;
                    index = index+2;
                end
            end
        end

        % Constraint 3 - V spatial consistency constraint
        disp('Generating inequality constraint 3');
        startIndex = 1;
        endIndex = startIndex + numWeightsPerPixel - 1;
        for t = 0:(T-2)
            for i = 1:M
                for j = 1:N
                    c = zeros(1, numVariables);
                    c(startIndex:endIndex) = beta;

                    % generate the mean computing matrices
                    startX = safeIndex(j, -windowSize, N);
                    endX = safeIndex(j, windowSize, N);
                    startY = safeIndex(i, -windowSize, M);
                    endY = safeIndex(i, windowSize, M);
                    xInterval = endX-startX+1;
                    yInterval = endY-startY+1;
                    nSummed = xInterval*yInterval;
                    startX = numWeightsPerPixel*(startX-1)+1;
                    endX = numWeightsPerPixel*endX;

                    f = zeros(M, numWeightsPerPixel*N);
                    f(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                        repmat(beta, [yInterval xInterval]);
                    f = reshape(f', [1 numWeightsPerFrame]);
                    g = zeros(1, numVariables);
                    g((t*numWeightsPerFrame+1):((t+1)*numWeightsPerFrame)) = f;

                    Z_index = numWeights + numPixels + (j-1)+(i-1)*N+t*M*N + 1;
                    Aineq(index+0, :) = c - g;
                    Aineq(index+0, Z_index) = -1;

                    Aineq(index+1, :) = g - c;
                    Aineq(index+1, Z_index) = -1;

                    startIndex = endIndex + 1;
                    endIndex = startIndex + numWeightsPerPixel - 1;
                    index = index+2;
                end
            end
        end

        % Constraint 4 - U momentum constraint
        disp('Generating inequality constraint 4');
        startIndex = 1;
        endIndex = startIndex + numWeightsPerPixel - 1;
        for t = 0:(T-3)
            for i = 1:M
                for j = 1:N
                    % generate the mean computing matrices
                    startX = safeIndex(j, -windowSize, N);
                    endX = safeIndex(j, windowSize, N);
                    startY = safeIndex(i, -windowSize, M);
                    endY = safeIndex(i, windowSize, M);
                    xInterval = endX-startX+1;
                    yInterval = endY-startY+1;
                    nSummed = xInterval*yInterval;
                    startX = numWeightsPerPixel*(startX-1)+1;
                    endX = numWeightsPerPixel*endX;

                    f = zeros(M, numWeightsPerPixel*N);
                    f(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                        repmat(alpha, [yInterval xInterval]);
                    f = reshape(f', [1 numWeightsPerFrame]);

                    g = zeros(1, numVariables);
                    g((t*numWeightsPerFrame+1):((t+1)*numWeightsPerFrame)) = f;
                    g(((t+1)*numWeightsPerFrame+1):((t+2)*numWeightsPerFrame)) = -f;

                    R_index = numWeights + numPixels + numPixels + ...
                        (j-1)+(i-1)*N+t*M*N + 1;
                    Aineq(index+0, :) = g;
                    Aineq(index+0, R_index) = -1;

                    Aineq(index+1, :) = -g;
                    Aineq(index+1, R_index) = -1;

                    startIndex = endIndex + 1;
                    endIndex = startIndex + numWeightsPerPixel - 1;
                    index = index+2;
                end
            end
        end

        % Constraint 5 - V momentum constraint
        disp('Generating inequality constraint 5');
        startIndex = 1;
        endIndex = startIndex + numWeightsPerPixel - 1;
        for t = 0:(T-3)
            for i = 1:M
                for j = 1:N
                    % generate the mean computing matrices
                    startX = safeIndex(j, -windowSize, N);
                    endX = safeIndex(j, windowSize, N);
                    startY = safeIndex(i, -windowSize, M);
                    endY = safeIndex(i, windowSize, M);
                    xInterval = endX-startX+1;
                    yInterval = endY-startY+1;
                    nSummed = xInterval*yInterval;
                    startX = numWeightsPerPixel*(startX-1)+1;
                    endX = numWeightsPerPixel*endX;

                    f = zeros(M, numWeightsPerPixel*N);
                    f(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                        repmat(beta, [yInterval xInterval]);
                    f = reshape(f', [1 numWeightsPerFrame]);
                    g = zeros(1, numVariables);
                    g((t*numWeightsPerFrame+1):((t+1)*numWeightsPerFrame)) = f;
                    g(((t+1)*numWeightsPerFrame+1):((t+2)*numWeightsPerFrame)) = -f;

                    S_index = numWeights + numPixels + numPixels + numMomentum + ...
                        (j-1)+(i-1)*N+t*M*N + 1;
                    Aineq(index+0, :) = g;
                    Aineq(index+0, S_index) = -1;

                    Aineq(index+1, :) = -g;
                    Aineq(index+1, S_index) = -1;

                    startIndex = endIndex + 1;
                    endIndex = startIndex + numWeightsPerPixel - 1;
                    index = index+2;
                end
            end
        end

        Aineq = sparse(Aineq);
        bineq = sparse(bineq);
    end
    
    if saveCSV
        save('weight_cache/Aeq.mat', 'Aeq');
        save('weight_cache/beq.mat', 'beq');
        save('weight_cache/Aineq.mat', 'Aineq');
        save('weight_cache/bineq.mat', 'bineq');
    end
    
    % Now, in theory, we should have the correct constraints to solve the QP
    W0 = zeros(1, numVariables);
    index = 1;
    for t = 1:(T-1)
        for i = 1:M
            for j = 1:N
                for b = 1:(2*windowSize+1)
                    for a = 1:(2*windowSize+1)
                        W0(index) = W{t}(i,j,a,b);
                        index = index+1;
                    end
                end
            end
        end
    end
    W0 = sparse(W0);

    if ~useL2Penalty
        disp('Solving linear program');
        [W_hat, fval, exitflag] = cplexlp(P, Aineq, bineq, Aeq, beq, ...
              [], [], W0);
    end
    
    % recover weights cell matrix from output vector
    newW = cell(1,T);
    index = 1;
    for t = 1:(T-1)
        newW{t} = zeros(M, N, win, win); 
        for i = 1:M
            for j = 1:N
                for b = 1:win
                    for a = 1:win
                        newW{t}(i, j, b, a) = W_hat(index);
                        index = index+1;
                    end
                end
            end
        end
    end
end


