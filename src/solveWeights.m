function [newW] = solveWeights(X, W, video, lambda, d)
%SOLVE_WEIGHTS Solve for the flow weights given a pixel labeling (using
% cplex)
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

    win = 2*d+1; % full window dimension
    [M, N] = size(X{1});
    T = video.T;
    numWeightsPerPixel = win^2;
    numPixelsPerFrame = M * N;
    numPixels = (T-1) * numPixelsPerFrame; % no weights for last frame
    numWeightsPerFrame = numWeightsPerPixel * numPixelsPerFrame;
    numWeights = numPixels * numWeightsPerPixel;
    numAux = (T-2) * numPixelsPerFrame;
    numVariables = numWeights + numPixels + numPixels + numAux + numAux;
    
    % T and F costs
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
                for b = (-d):d
                    for a = (-d):d
                        i_tplus = safeIndex(i, a, M);
                        j_tplus = safeIndex(j, b, N);
                       
                        P(index) = abs(X_t(i,j) - X_tplus(i_tplus,j_tplus));
%                         Q(index) = abs(I_t(i,j,1) - I_tplus(i_tplus,j_tplus,1)) + ...
%                                 abs(I_t(i,j,2) - I_tplus(i_tplus,j_tplus,2)) + ...
%                                 abs(I_t(i,j,3) - I_tplus(i_tplus,j_tplus,3));
                        Q(index) = abs(I_t(i,j) - I_tplus(i_tplus,j_tplus));
                        
                        index = index+1;
                    end
                end
            end
        end
    end
    P = sparse(P);
    Q = sparse(Q);
    
    % C cost - create vector to sum auxiliary flow difference variables
    C = zeros(numVariables,1);
    C((numWeights+1):(numWeights+2*numPixels)) = ones(2*numPixels,1);
    C = sparse(C);
    
    % M cost - create H matrix for quadratic momentum penalty
    H = zeros(numVariables, numVariables);
    if numAux > 0
        H(1:numAux,(numWeights+2*numPixels+1):numVariables) = ...
            [eye(numAux) eye(numAux)];
        H((numWeights+2*numPixels+1):numVariables,1:numAux) = ...
            [eye(numAux); eye(numAux)];
    end
    H = lambda(6)*H;
    
    P = lambda(3)*P + lambda(4)*Q + lambda(5)*C; % final linear cost term
    
    % Equality constraints - ensure weights sum to 1
    Aeq = zeros(numPixels, numVariables);
    startIndex = 1;
    endIndex = startIndex + numWeightsPerPixel -1;
    for i = 1:numPixels
       Aeq(i,startIndex:endIndex) = ones(1, numWeightsPerPixel);
       startIndex = endIndex + 1;
       endIndex = startIndex + numWeightsPerPixel - 1;
    end
    beq = ones(numPixels, 1);
    
    % Inequality constraints
    % 1. Weights greater than 0
    % 2. Y auxiliary variables are greater that +/- the U flow continuity
    % 3. Z auxiliary variables are greater that +/- the V flow continuity
    % 4. R auxiliary variables are greater that +/- the U flow momentum
    % 5. S auxiliary variables are greater that +/- the V flow momentum
    numInequalityConstraints = numWeights + 2*numPixels + 2*numPixels + ...
        2*numAux + 2*numAux;
    Aineq = sparse(numInequalityConstraints, numVariables);
    bineq = sparse(numInequalityConstraints, 1);
    
    % Constraint 1
    Aineq(1:numWeights,1:numWeights) = -eye(numWeights);
    index = numWeights+1;
    
    % helper vectors for computing U, V
    alpha = repmat(linspace(-d, d, win), [1, win]);
    beta = reshape(repmat(linspace(-d, d, win), [win, 1]), [1 win^2]);
    
    % Constraint 2 - U spatial consistency constraint
    startIndex = 1;
    endIndex = startIndex + numWeightsPerPixel - 1;
    for t = 0:(T-2)
        for i = 1:M
            for j = 1:N
                c = sparse(1, numVariables);
                c(startIndex:endIndex) = alpha;

                % generate the mean computing matrices
                startX = safeIndex(j, -d, N);
                endX = safeIndex(j, d, N);
                startY = safeIndex(i, -d, M);
                endY = safeIndex(i, d, M);
                xInterval = endX-startX+1;
                yInterval = endY-startY+1;
                nSummed = xInterval*yInterval;
                startX = numWeightsPerPixel*(startX-1)+1;
                endX = numWeightsPerPixel*endX;
                
                f = sparse(M, numWeightsPerPixel*N);
                f(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    repmat(alpha, [yInterval xInterval]);
                f = reshape(f', [1 numWeightsPerFrame]);
                g = sparse(1, numVariables);
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
    startIndex = 1;
    endIndex = startIndex + numWeightsPerPixel - 1;
    for t = 0:(T-2)
        for i = 1:M
            for j = 1:N
                c = sparse(1, numVariables);
                c(startIndex:endIndex) = beta;

                % generate the mean computing matrices
                startX = safeIndex(j, -d, N);
                endX = safeIndex(j, d, N);
                startY = safeIndex(i, -d, M);
                endY = safeIndex(i, d, M);
                xInterval = endX-startX+1;
                yInterval = endY-startY+1;
                nSummed = xInterval*yInterval;
                startX = numWeightsPerPixel*(startX-1)+1;
                endX = numWeightsPerPixel*endX;
                
                f = sparse(M, numWeightsPerPixel*N);
                f(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    repmat(beta, [yInterval xInterval]);
                f = reshape(f', [1 numWeightsPerFrame]);
                g = sparse(1, numVariables);
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
    for t = 0:(T-2)
        for i = 1:M
            for j = 1:N
                for b = (-d):d
                    for a = (-d):d
  
                        i_tplus = max(min(i+b, M), 1);
                        j_tplus = max(min(j+a, N), 1);
                        startIndex = (j_tplus-1)*numWeightsPerPixel + ...
                            (i_tplus-1)*numWeightsPerPixel*N + ...
                            (t+1)*numWeightsPerPixel*N*M + 1;
                        endIndex = startIndex + numWeightsPerPixel - 1;
                        
                        c = sparse(1, numVariables);
                        c(startIndex:endIndex) = alpha;
                        
                        R_index = numWeights + 2*numPixels + ...
                            (a+d) + (b+d)*win + ...
                            (j-1)*numWeightsPerPixel + ...
                            (i-1)*numWeightsPerPixel*N + ...
                            t*numWeightsPerPixel*N*M + 1;
                        Aineq(index+0, :) = -c;
                        Aineq(index+0, R_index) = -1;

                        Aineq(index+1, :) = c;
                        Aineq(index+1, R_index) = -1;

                        bineq(index+0) = -a;
                        bineq(index+1) = a;
                        
                        index = index+2;
                    end
                end
            end
        end
    end
    
    % Constraint 5 - V momentum constraint
    for t = 0:(T-2)
        for i = 1:M
            for j = 1:N
                for b = (-d):d
                    for a = (-d):d
  
                        i_tplus = max(min(i+b, M), 1);
                        j_tplus = max(min(j+a, N), 1);
                        startIndex = (j_tplus-1)*numWeightsPerPixel + ...
                            (i_tplus-1)*numWeightsPerPixel*N + ...
                            (t+1)*numWeightsPerPixel*N*M + 1;
                        endIndex = startIndex + numWeightsPerPixel - 1;
                        
                        c = sparse(1, numVariables);
                        c(startIndex:endIndex) = beta;
                        
                        S_index = numWeights + 2*numPixels + numAux + ...
                            (a+d) + (b+d)*win + ...
                            (j-1)*numWeightsPerPixel + ...
                            (i-1)*numWeightsPerPixel*N + ...
                            t*numWeightsPerPixel*N*M + 1;
                        Aineq(index+0, :) = -c;
                        Aineq(index+0, S_index) = -1;

                        Aineq(index+1, :) = c;
                        Aineq(index+1, S_index) = -1;

                        bineq(index+0) = -b;
                        bineq(index+1) = b;
                        
                        index = index+2;
                    end
                end
            end
        end
    end
    
    % Now, in theory, we shoudl have the correct constraints to solve the QP
    W0 = sparse(1, numVariables);
    index = 1;
    for t = 1:(T-1)
        for i = 1:M
            for j = 1:N
                for b = 1:(2*d+1)
                    for a = 1:(2*d+1)
                        W0(index) = W{t}(i,j,a,b);
                        index = index+1;
                    end
                end
            end
        end
    end

%     [W_hat, fval, exitflag] = cplexlp(P, Aineq, bineq, Aeq, beq, ...
%         [], [], W0);
     [W_hat, fval, exitflag] = cplexqp(H, P, Aineq, bineq, Aeq, beq, ...
         [], [], W0);
    
    % recover weights cell matrix from output vector
    newW = cell(1,video.T);
    startIndex = 1;
    endIndex = startIndex + numWeightsPerFrame;
    for t = 1:(T-1)
        newW{t} = W_hat(startIndex:endIndex);
        startIndex = endIndex + 1;
        endIndex = startIndex + numWeightsPerFrame;
    end
    newW{T} = W{T};
    
    for i = 1:numPixels
       if W_hat((i-1)*numWeightsPerPixel+1) ~= 0
           fprintf('Error!\n');
       end
    end
end

