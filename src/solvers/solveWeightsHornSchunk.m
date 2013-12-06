function [newW] = solveWeightsHornSchunk(X, U, V, video, T, lambda, d)
%SOLVE_WEIGHT_AVG_MOMENTUM Solve for the flow weights given a pixel labeling (using
% cplex) with average momentum constraint
% The variable we end up optimizing is W_hat = [W Y Z R S]', where
% Y puts the U spatial constraints into epigraph form,
% Z puts the V spatial constraints into epigraph form,
% R puts the U momentum constraints into epigraph form, and
% S puts the V momentum constraints into epigraph form
%   INPUTS -
%   X : tentative labels of pixels in the video sequence
%   U : initial U components of flow
%   V : initial V components of flow
%   video : struct containing the video data
%   lambda : vector of constant weights in the objective function
%   d : size of window to average flows over
%
% OUTPUTS - 
%   returns the weights found by cplex
%   example : generateSyntheticDataMovement(seedImg,2+2*sin(1:20),2+2*cos(1:20),20,0.2);

    [M, N] = size(X{1});
    c = size(video.I, 3);
    numPixelsPerFrame = M * N;
    numFlowsPerFrame = 2*numPixelsPerFrame;
    numFlows= (T-1) * numFlowsPerFrame; % no weights for last frame
    numMomentum = 2 * (T-2) * numPixelsPerFrame;
    numVariables = 4*numFlows + c*numFlows + 2*numMomentum;

    disp('Creating linear costs');
    
    % T cost - create vector to sum auxiliary first order label constancy variables
    T = zeros(numVariables,1);
    T((numFlows+1):(2*numFlows)) = ones(numFlows,1);
    T = sparse(T);
    
    % F cost - create vector to sum auxiliary first order image intensity constancy variables
    F = zeros(numVariables,1);
    F((2*numFlows+1):(3*numFlows)) = ones(numFlows,1); % grayscale only for now
    F = sparse(F);
    
    % C cost - create vector to sum auxiliary flow difference variables
    C = zeros(numVariables,1);
    C((3*numFlows+1):(5*numFlows)) = ones(2*numFlows,1);
    C = sparse(C);
    
    % M cost - create vector to sum auxiliary momentum variables
    H = zeros(numVariables,1);
    H((5*numFlows+1):numVariables) = ones(2*numMomentum,1);
    H = sparse(H);

    T = lambda(3)*T + lambda(4)*F + lambda(5)*C + lambda(6)*H; % final linear cost term
    
    % Inequality constraints
    % 1. P auxiliary variables are +/- the spatial similarity (for the
    % label image)
    % 2. Q auxiliary variables are +/- the spatial similarity (for the
    % intensity / color image)
    % 2. Y auxiliary variables are greater than +/- the U flow continuity
    % 3. Z auxiliary variables are greater than +/- the V flow continuity
    % 4. R auxiliary variables are greater than +/- the U avg flow diff
    % between frames
    % 5. S auxiliary variables are greater than +/- the V avg flow diff
    % between frames
    disp('Generating inequality constraints');
    numInequalityConstraints = 2*(6*numPixels);
    Aineq = zeros(numInequalityConstraints, numVariables);
    bineq = zeros(numInequalityConstraints, 1);
    
    % helper functions for indexing into the arrays
    linearIndex = @(i,j,t) ((j-1) + 2*(i-1)*N + 2*(t-1)*M*N + 1);
    linearUIndex = @(i,j,t) ((j-1) + 2*(i-1)*N + 2*(t-1)*M*N + 1);
    linearVIndex = @(i,j,t) ((j-1+N) + 2*(i-1)*N + 2*(t-1)*M*N + 1);
       
    % Constraint 1 - P label similarity constraint
    disp('Generating inequality constraint 1');
    index = 1;  
    for t = 1:(T-1)
        [labelGradientX, labelGradientY] = imgradientxy(X{t});
        for i = 1:M
            for j = 1:N
                pixelLinearIndex = linearIndex(i,j,t);
                U_index = linearUIndex(i,j,t);
                V_index = linearVIndex(i,j,t);
                P_index = numFlows+pixelLinearIndex;
                labelGradientT = X{t+1}(i,j) - X{t}(i,j);
                
                % add positive inequality constraint
                Aineq(index+0, U_index) = labelGradientX(i,j);
                Aineq(index+0, V_index) = labelGradientY(i,j);
                Aineq(index+0, P_index) = -1;
                bineq(index+0) = -labelGradientT;
                
                % add negative inequality constraint
                Aineq(index+1, U_index) = -labelGradientX(i,j);
                Aineq(index+1, V_index) = -labelGradientY(i,j);
                Aineq(index+1, P_index) = -1;
                bineq(index+1) = labelGradientT;

                index = index+2;
            end
        end
    end
    
    % Constraint 2 - Q image intensity similarity constraint
    disp('Generating inequality constraint 2');
    index = 1;  
    for t = 1:(T-1)
        [imGradientX, imGradientY] = imgradientxy(video.I{t});
        for i = 1:M
            for j = 1:N
                pixelLinearIndex = linearIndex(i,j,t);
                U_index = linearUIndex(i,j,t);
                V_index = linearVIndex(i,j,t);
                Q_index = 2*numFlows + pixelLinearIndex;
                imGradientT = video.I{t+1}(i,j) - video.I{t}(i,j);
                
                % add positive inequality constraint
                Aineq(index+0, U_index) = imGradientX(i,j);
                Aineq(index+0, V_index) = imGradientY(i,j);
                Aineq(index+0, Q_index) = -1;
                bineq(index+0) = -imGradientT;
                
                % add negative inequality constraint
                Aineq(index+1, U_index) = -imGradientX(i,j);
                Aineq(index+1, V_index) = -imGradientY(i,j);
                Aineq(index+1, Q_index) = -1;
                bineq(index+1) = imGradientT;

                index = index+2;
            end
        end
    end
    
    % Constraint 3 - Y - U spatial consistency constraint
    disp('Generating inequality constraint 3');
    startIndex = 1;
    endIndex = startIndex + numFlows - 1;
    for t = 1:(T-1)
        for i = 1:M
            for j = 1:N
                % generate the mean computing matrices
                startX = safeIndex(j, -d, N);
                endX = safeIndex(j, d, N);
                startY = safeIndex(i, -d, M);
                endY = safeIndex(i, d, M);
                xInterval = endX-startX+1;
                yInterval = endY-startY+1;
                nSummed = xInterval*yInterval;
                
                avgFilter = zeros(M, N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    ones(yInterval, xInterval);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);
                g = zeros(1, numVariables);
                g((t*numFlowsPerFrame+1):((t+1)*numFlowsPerFrame-numPixelsPerFrame)) = ...
                    avgFilter;
                
                pixelLinearIndex = linearIndex(i,j,t);
                U_index = linearUIndex(i,j,t);
                Y_index = 3*numFlows + pixelLinearIndex;
                
                Aineq(index+0, :) = -g;
                Aineq(index+0, U_index) = 1;
                Aineq(index+0, Y_index) = -1;
                
                Aineq(index+1, :) = g;
                Aineq(index+1, U_index) = -1;
                Aineq(index+1, Y_index) = -1;

                index = index+2;
            end
        end
    end
    
    % Constraint 4 - Z - V spatial consistency constraint
    disp('Generating inequality constraint 4');
    startIndex = 1;
    endIndex = startIndex + numFlows - 1;
    for t = 1:(T-1)
        for i = 1:M
            for j = 1:N
                % generate the mean computing matrices
                startX = safeIndex(j, -d, N);
                endX = safeIndex(j, d, N);
                startY = safeIndex(i, -d, M);
                endY = safeIndex(i, d, M);
                xInterval = endX-startX+1;
                yInterval = endY-startY+1;
                nSummed = xInterval*yInterval;
                
                avgFilter = zeros(M, N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    ones(yInterval, xInterval);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);
                g = zeros(1, numVariables);
                g((t*numFlowsPerFrame+numPixelsPerFrame+1):((t+1)*numFlowsPerFrame)) = ...
                    avgFilter;
                
                pixelLinearIndex = linearIndex(i,j,t);
                V_index = linearVIndex(i,j,t);
                Z_index = 4*numFlows + pixelLinearIndex;
                
                Aineq(index+0, :) = -g;
                Aineq(index+0, V_index) = 1;
                Aineq(index+0, Z_index) = -1;
                
                Aineq(index+1, :) = g;
                Aineq(index+1, V_index) = -1;
                Aineq(index+1, Z_index) = -1;

                index = index+2;
            end
        end
    end
    
    % Constraint 5 - R - U momentum constraint
    disp('Generating inequality constraint 4');
    startIndex = 1;
    endIndex = startIndex + numWeightsPerPixel - 1;
    for t = 1:(T-2)
        for i = 1:M
            for j = 1:N
                % generate the mean computing matrices
                startX = safeIndex(j, -d, N);
                endX = safeIndex(j, d, N);
                startY = safeIndex(i, -d, M);
                endY = safeIndex(i, d, M);
                xInterval = endX-startX+1;
                yInterval = endY-startY+1;
                nSummed = xInterval*yInterval;
                
                avgFilter = zeros(M, numWeightsPerPixel*N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    repmat(alpha, [yInterval xInterval]);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);

                g = zeros(1, numVariables);
                g((t*numFlowsPerFrame+1):((t+1)*numFlowsPerFrame-numPixelsPerFrame)) = ... 
                    avgFilter;
                g(((t+1)*numFlowsPerFrame+1):((t+2)*numFlowsPerFrame-numPixelsPerFrame)) = ...
                    -avgFilter;

                pixelLinearIndex = linearIndex(i,j,t);
                R_index = 5*numPixels + pixelLinearIndex;
                
                Aineq(index+0, :) = g;
                Aineq(index+0, R_index) = -1;
                
                Aineq(index+1, :) = -g;
                Aineq(index+1, R_index) = -1;

                index = index+2;
            end
        end
    end
    
    % Constraint 6 - S - V momentum constraint
    disp('Generating inequality constraint 6');
    startIndex = 1;
    endIndex = startIndex + numWeightsPerPixel - 1;
    for t = 1:(T-2)
        for i = 1:M
            for j = 1:N
                % generate the mean computing matrices
                startX = safeIndex(j, -d, N);
                endX = safeIndex(j, d, N);
                startY = safeIndex(i, -d, M);
                endY = safeIndex(i, d, M);
                xInterval = endX-startX+1;
                yInterval = endY-startY+1;
                nSummed = xInterval*yInterval;
                
                avgFilter = zeros(M, numWeightsPerPixel*N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    repmat(alpha, [yInterval xInterval]);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);

                g = zeros(1, numVariables);
                g((t*numFlowsPerFrame+numPixelsPerFrame+1):((t+1)*numFlowsPerFrame)) = ... 
                    avgFilter;
                g(((t+1)*numFlowsPerFrame+numPixelsPerFrame+1):((t+2)*numFlowsPerFrame)) = ...
                    -avgFilter;

                pixelLinearIndex = linearIndex(i,j,t);
                R_index = 6*numPixels + pixelLinearIndex;
                
                Aineq(index+0, :) = g;
                Aineq(index+0, R_index) = -1;
                
                Aineq(index+1, :) = -g;
                Aineq(index+1, R_index) = -1;

                index = index+2;
            end
        end
    end  
    Aineq = sparse(Aineq);
    bineq = sparse(bineq);

    disp('Solving linear program');
    [newUV, fval, exitflag] = cplexlp(T, Aineq, bineq, Aeq, beq, ...
          [], [], W0);
end


