function [newU, newV] = solveWeightsHornSchunk(X, U, V, video, T, lambda, d, useL2Penalty, debug)
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
    numPixels = (T-1) * numPixelsPerFrame;
    numFlows= (T-1) * numFlowsPerFrame; % no weights for last frame
    numMomentum = (T-2) * numPixelsPerFrame;
    numVariables = numFlows + numPixels + c*numPixels + numFlows + ...
        2*numMomentum;
    numVariableTypes = 8;

    % create double versions of images
    I = cell(T);
    for t = 1:T
       I{t} = double(video.I{t}); 
    end
    
    disp('Creating linear costs');
    
    % Compute all of the start indices
    startIndices = zeros(numVariableTypes,1);
    startIndices(1) = 1;                            % U
    startIndices(2) = startIndices(1) + numPixels;  % V
    startIndices(3) = startIndices(2) + numPixels;  % P
    startIndices(4) = startIndices(3) + numPixels;  % Q
    startIndices(5) = startIndices(4) + numPixels;  % Y
    startIndices(6) = startIndices(5) + numPixels;  % Z
    startIndices(7) = startIndices(6) + numPixels;  % R
    startIndices(8) = startIndices(7) + numMomentum;% S
    
    % Create cost vector to sum up the auxiliary variables
    L = zeros(numVariables,1);
    L(startIndices(3):startIndices(4)-1) = lambda(3) * ones(numPixels,1); % P summation
    L(startIndices(4):startIndices(5)-1) = lambda(4) * ones(numPixels,1); % Q summation (grayscale only for now)
    L(startIndices(5):startIndices(6)-1) = lambda(5) * ones(numPixels,1);  % Y summation
    L(startIndices(6):startIndices(7)-1) = lambda(5) * ones(numPixels,1);  % Z summation
    L(startIndices(7):startIndices(8)-1) = lambda(6) * ones(numMomentum,1); % R summation
    L(startIndices(8):numVariables) = lambda(6) * ones(numMomentum,1); % S summation
    L = sparse(L);
    
    % create diagonal quadratic penalty if specified
    H = [];
    if useL2Penalty
        H = lambda(7)*eye(numVariables);
    end
   
    % Inequality constraints
    % 1. P auxiliary variables are +/- the spatial similarity (for the
    % label image) - Horn-Shunk method
    % 2. Q auxiliary variables are +/- the spatial similarity (for the
    % intensity / color image) - Horn-Shunk method
    % 2. Y auxiliary variables are greater than +/- the U flow continuity
    % 3. Z auxiliary variables are greater than +/- the V flow continuity
    % 4. R auxiliary variables are greater than +/- the U avg flow diff
    % between frames
    % 5. S auxiliary variables are greater than +/- the V avg flow diff
    % between frames
    % 6. Lower bound on UV vars
    % 7. Upper bound on UV vars
    disp('Generating inequality constraints');
    numInequalityConstraints = 2*(2*numPixels) + 2*(4*numFlows) + 2*numFlows;
    Aineq = zeros(numInequalityConstraints, numVariables);
    bineq = zeros(numInequalityConstraints, 1);
    
    % helper functions for indexing into the arrays
    linearIndex = @(i,j,t) ((j-1) + (i-1)*N + (t-1)*M*N + 1);
    linearUIndex = @(i,j,t) (startIndices(1) + linearIndex(i,j,t) - 1);
    linearVIndex = @(i,j,t) (startIndices(2) + linearIndex(i,j,t) - 1);
    linearPIndex = @(i,j,t) (startIndices(3) + linearIndex(i,j,t) - 1);
    linearQIndex = @(i,j,t) (startIndices(4) + linearIndex(i,j,t) - 1);
    linearYIndex = @(i,j,t) (startIndices(5) + linearIndex(i,j,t) - 1);
    linearZIndex = @(i,j,t) (startIndices(6) + linearIndex(i,j,t) - 1);
    linearRIndex = @(i,j,t) (startIndices(7) + linearIndex(i,j,t) - 1);
    linearSIndex = @(i,j,t) (startIndices(8) + linearIndex(i,j,t) - 1);
       
    % Constraint 1 - P label similarity constraint
    disp('Generating inequality constraint 1');
    index = 1; % keeps track of the current constraint for ease of coding
    for t = 1:(T-1)
        [labelGradientX, labelGradientY] = imgradientxy(X{t});
        for i = 1:M
            for j = 1:N
                U_index = linearUIndex(i,j,t);
                V_index = linearVIndex(i,j,t);
                P_index = linearPIndex(i,j,t);
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
    for t = 1:(T-1)
        [imGradientX, imGradientY] = imgradientxy(I{t});
        for i = 1:M
            for j = 1:N
                U_index = linearUIndex(i,j,t);
                V_index = linearVIndex(i,j,t);
                Q_index = linearQIndex(i,j,t);
                imGradientT = I{t+1}(i,j) - I{t}(i,j);
                
                % add positive inequality constraint
                g = 1.0; % constant to normalize image gradients (NOT USED)
                Aineq(index+0, U_index) = g*imGradientX(i,j);
                Aineq(index+0, V_index) = g*imGradientY(i,j);
                Aineq(index+0, Q_index) = -1;
                bineq(index+0) = -g*imGradientT;
                
                % add negative inequality constraint
                Aineq(index+1, U_index) = -g*imGradientX(i,j);
                Aineq(index+1, V_index) = -g*imGradientY(i,j);
                Aineq(index+1, Q_index) = -1;
                bineq(index+1) = g*imGradientT;

                index = index+2;
            end
        end
    end
    
    % Constraint 3 - Y - U spatial consistency constraint
    disp('Generating inequality constraint 3');
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
                
                % build an averaging filter by first putting the summation
                % window into the location in the image and then reshaping
                % to match the size of the variable vector 
                avgFilter = zeros(M, N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    ones(yInterval, xInterval);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);
                
                avgFilterAllVariables = zeros(1, numVariables);
                startIndex = linearUIndex(1,1,t);
                endIndex = linearUIndex(1,1,t+1)-1;
                avgFilterAllVariables(startIndex:endIndex) = avgFilter;
                
                U_index = linearUIndex(i,j,t);
                Y_index = linearYIndex(i,j,t);
                
                % place vectors into inequality matrices
                Aineq(index+0, :) = -avgFilterAllVariables;
                Aineq(index+0, U_index) = Aineq(index+0, U_index) + 1;
                Aineq(index+0, Y_index) = -1;
                
                Aineq(index+1, :) = avgFilterAllVariables;
                Aineq(index+1, U_index) = Aineq(index+1, U_index) - 1;
                Aineq(index+1, Y_index) = -1;

                index = index+2;
            end
        end
    end
    
    % Constraint 4 - Z - V spatial consistency constraint
    disp('Generating inequality constraint 4');
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
                
                % build an averaging filter by first putting the summation
                % window into the location in the image and then reshaping
                % to match the size of the variable vector 
                avgFilter = zeros(M, N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    ones(yInterval, xInterval);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);
                
                avgFilterAllVariables = zeros(1, numVariables);
                startIndex = linearVIndex(1,1,t);
                endIndex = linearVIndex(1,1,t+1)-1;
                avgFilterAllVariables(startIndex:endIndex) = avgFilter;
                
                V_index = linearVIndex(i,j,t);
                Z_index = linearZIndex(i,j,t);
               
                % place vectors into inequality matrices
                Aineq(index+0, :) = -avgFilterAllVariables;
                Aineq(index+0, V_index) = Aineq(index+0, V_index) + 1;
                Aineq(index+0, Z_index) = -1;
                
                Aineq(index+1, :) = avgFilterAllVariables;
                Aineq(index+1, V_index) = Aineq(index+1, V_index) - 1;
                Aineq(index+1, Z_index) = -1;

                index = index+2;
            end
        end
    end
    
    % Constraint 5 - R - U momentum constraint
    disp('Generating inequality constraint 5');
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
                
                % build an averaging filter by first putting the summation
                % window into the location in the image and then reshaping
                % to match the size of the variable vectors
                avgFilter = zeros(M, N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    ones(yInterval, xInterval);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);

                % put the filter in the variable range for times t and t+1
                % since we want to compute the temporal derivatives
                avgFilterAllVariables = zeros(1, numVariables);
                startIndex = linearUIndex(1,1,t);
                endIndex = linearUIndex(1,1,t+1)-1;
                avgFilterAllVariables(startIndex:endIndex) = avgFilter;
                startIndex = linearUIndex(1,1,t+1);
                endIndex = linearUIndex(1,1,t+2)-1;
                avgFilterAllVariables(startIndex:endIndex) = -avgFilter;
                
                R_index = linearRIndex(i,j,t);               
                
                % update inequalities
                Aineq(index+0, :) = avgFilterAllVariables;
                Aineq(index+0, R_index) = -1;
                
                Aineq(index+1, :) = -avgFilterAllVariables;
                Aineq(index+1, R_index) = -1;

                index = index+2;
            end
        end
    end
    
    % Constraint 6 - S - V momentum constraint
    disp('Generating inequality constraint 6');
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
                
                % build an averaging filter by first putting the summation
                % window into the location in the image and then reshaping
                % to match the size of the variable vectors
                avgFilter = zeros(M, N);
                avgFilter(startY:endY, startX:endX) = (1.0 / nSummed) * ...
                    ones(yInterval, xInterval);
                avgFilter = reshape(avgFilter', [1 numPixelsPerFrame]);

                % put the filter in the variable range for times t and t+1
                % since we want to compute the temporal derivatives
                avgFilterAllVariables = zeros(1, numVariables);
                startIndex = linearVIndex(1,1,t);
                endIndex = linearVIndex(1,1,t+1)-1;
                avgFilterAllVariables(startIndex:endIndex) = avgFilter;
                startIndex = linearVIndex(1,1,t+1);
                endIndex = linearVIndex(1,1,t+2)-1;
                avgFilterAllVariables(startIndex:endIndex) = -avgFilter;

                S_index = linearSIndex(i,j,t);               
                
                Aineq(index+0, :) = avgFilterAllVariables;
                Aineq(index+0, S_index) = -1;
                
                Aineq(index+1, :) = -avgFilterAllVariables;
                Aineq(index+1, S_index) = -1;

                index = index+2;
            end
        end
    end 
    
    % Constraints 7,8 - Lower/upper bounds on UV
    disp('Generating inequality constraint 7');
    startIndex = index;
    endIndex = startIndex + numFlows - 1;
    Aineq(startIndex:endIndex, 1:numFlows) = eye(numFlows);
    bineq(startIndex:endIndex) = d*ones(numFlows,1);
    index = index + numFlows + 1;
    
    disp('Generating inequality constraint 8');
    startIndex = index;
    endIndex = startIndex + numFlows - 1;
    Aineq(startIndex:endIndex,1:numFlows) = -eye(numFlows);
    bineq(startIndex:endIndex) = d*ones(numFlows,1);
    index = index + numFlows + 1;
    
    Aineq = sparse(Aineq);
    bineq = sparse(bineq);

    % Create the estimated flow vector
    UV_init = zeros(numVariables, 1);
    for t = 1:(T-2)
        for i = 1:M
            for j = 1:N
                U_index = linearUIndex(i,j,t);
                V_index = linearVIndex(i,j,t);
                UV_init(U_index) = U{t}(i,j);
                UV_init(V_index) = V{t}(i,j);
            end
        end
    end
    
    if useL2Penalty
        disp('Solving quadratic program');
        [UV_new, fval, exitflag, output] = cplexqp(H, L, Aineq, bineq, [], [], ...
          [], [], UV_init);
    else
        disp('Solving linear program');
        [UV_new, fval, exitflag, output] = cplexlp(L, Aineq, bineq, [], [], ...
          [], [], UV_init);
    end
    
    if debug
        % query pixel index for debugging
        i = 6;
        j = 6;
        t = 1;
        
        % pull variable values for spcified index
        ui = linearUIndex(i,j,t);
        vi = linearVIndex(i,j,t);
        pi = linearPIndex(i,j,t);
        qi = linearQIndex(i,j,t);
        yi = linearYIndex(i,j,t);
        u = UV_new(ui);
        v = UV_new(vi);
        p = UV_new(pi);
        q = UV_new(qi);
        y = UV_new(yi);
        Xt = X{t+1}(i,j) - X{t}(i,j);
        [Xxp, Xyp] = imgradientxy(X{t});
        Xx = Xxp(i,j);
        Xy = Xyp(i,j);
        It = I{t+1}(i,j) - I{t}(i,j);
        [Ixp, Iyp] = imgradientxy(I{t});
        Ix = Ixp(i,j);
        Iy = Iyp(i,j);
        test = 1;
    end
          
    % recover flows cell matrix from output vector
    newU = cell(1,T);
    newV = cell(1,T);
    index = 1;
    for t = 1:(T-1)
        newU{t} = zeros(M, N); 
        newV{t} = zeros(M, N); 
        for i = 1:M
            for j = 1:N
                % here we assume all t of U come before all t of V,
                % but above its U,V alternating!
                newU{t}(i, j) = UV_new(startIndices(1)+index-1);
                newV{t}(i, j) = UV_new(startIndices(2)+index-1);
                index = index+1;
            end
        end
    end
    newU{T} = U{T};
    newV{T} = V{T};
end


