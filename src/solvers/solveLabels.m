function [newX] = solveLabels(X, U, V, video, T, ...
    lambda, spatialWinSize, temporalWinSize, useL2Penalty, ...
    loadCSV, saveCSV, debug)
%SOLVE_WEIGHT_AVG_MOMENTUM Solve for the flow weights given a pixel labeling (using
% cplex) with average momentum constraint
% The variable we end up optimizing is W_hat = [X Y Z]', where
% Y puts the X spatial constraints into epigraph form,
% Z puts the X temporal constraints into epigraph form,
% Yes, we're still missing the size constraint because it's less important
%   INPUTS -
%   X : tentative labels of pixels in the video sequence
%   U : initial U components of flow
%   V : initial V components of flow
%   video : struct containing the video data
%   lambda : vector of constant weights in the objective function
%   windowSize : size of window to average flows over
%   useL2Penalty: whether or not to use an L2 penalty
%   loadCSV: whether or not to load the constraint matrix from a CSV
%   saveCSV: whether or not to save the constraint matrix to a CSV
%   debug: flag to stop program in an if statement to examine specific
%   variables
%
% OUTPUTS - 
%   returns the weights found by cplex
%   example : generateSyntheticDataMovement(seedImg,2+2*sin(1:20),2+2*cos(1:20),20,0.2);

    [M, N] = size(X{1});
    c = size(video.I, 3);
    spatialWinDim = 2*spatialWinSize+1;
    temporalWinDim = 2*temporalWinSize+1;
    
    numPixelsPerFrame = M * N;
    numFlowsPerFrame = 2*numPixelsPerFrame;
    numWeightsPerFrame =  temporalWinDim^2 * numPixelsPerFrame;
    numLabels = T * numPixelsPerFrame;
    numFlows= (T-1) * numFlowsPerFrame;
    numSpatial = spatialWinDim^2 * numLabels;
    numTemporal = (T-1) * numWeightsPerFrame;
    numVariables = numLabels + numSpatial + numTemporal;
    numVariableTypes = 3;

    % create double versions of images
    I = cell(T);
    for t = 1:T
       I{t} = double(video.I{t}); 
    end
    
    % get weights
    W_vec = zeros(numTemporal,1);
    startIndex = 1;
    endIndex = startIndex + temporalWinDim^2 - 1;
    temp = zeros(temporalWinDim, temporalWinDim);
    for t = 1:(T-1)
        W{t} = uv_to_weights(U{t}, V{t}, temporalWinSize);
        for i = 1:M
            for j = 1:N
                temp(:,:) = W{t}(i,j,:,:);
                W_vec(startIndex:endIndex) = reshape(temp',...
                    [temporalWinDim^2, 1]);
                startIndex = endIndex+1;
                endIndex = startIndex + temporalWinDim^2 - 1;
            end
        end
    end
    disp('Creating linear costs');
    
    % Compute all of the start indices
    startIndices = zeros(numVariableTypes,1);
    startIndices(1) = 1;                            % X
    startIndices(2) = startIndices(1) + numLabels;  % Y
    startIndices(3) = startIndices(2) + numSpatial; % Z
    
    % TODO: Create the size change variable n stuff
    
    % normalize lambdas
    lambda(1) = double(lambda(1)) / double(M*N*T);
    lambda(2) = double(lambda(2)) / double(M*N*T*spatialWinDim^2);
    lambda(3) = double(lambda(3)) / double(M*N*(T-1));
    
    % Get cost acc. to appearance model for each frame t in T
    C = randomForestBasedCost(video);

    % if C is a cell make it a C(T,M,N)
    if iscell (C)
        for t = 1:T
            Ctemp(t,:,:) = C{t};
        end
    end
    C = Ctemp;
    clear Ctemp

    % rewrite Cost matrix C(T,M,N) into C(T, M*N)
    appearanceCost = reshape(C,T, M*N);
    appearanceCost = appearanceCost';
    appearanceCost = reshape(appearanceCost, [numLabels, 1]);
    
    % Create cost vector to sum up the auxiliary variables
    if loadCSV
        L = load('label_cache/L.mat', 'L');
        L = L.L;

        % create diagonal quadratic penalty if specified
        if useL2Penalty
            H = load('label_cache/H.mat', 'H');
            H = H.H;
        end
    else
        L = zeros(numVariables,1);
        L(startIndices(1):startIndices(2)-1) = appearanceCost; % X summation
        L(startIndices(2):startIndices(3)-1) = ones(numSpatial,1); % Y summation
        L(startIndices(3):numVariables) = W_vec;     % Z summation    

        % create diagonal quadratic penalty if specified
        H = [];
        if useL2Penalty
            H = lambda(7)*eye(numVariables);
        end
    end
    C = zeros(numVariables,1);   
    C(startIndices(1):startIndices(2)-1) = lambda(1) * L(startIndices(1):startIndices(2)-1); % X summation
    C(startIndices(2):startIndices(3)-1) = lambda(2) * L(startIndices(2):startIndices(3)-1); % Y summation
    C(startIndices(3):numVariables) = lambda(3) * W_vec;     % Z summation    
    C = sparse(C);

   
    % Inequality constraints
    % 1. Y auxiliary variables are +/- the spatial similarity (for the
    % label image)
    % 2. Z auxiliary variables are +/- the temporal similarity (for the
    % label image)
    disp('Generating inequality constraints');
    numInequalityConstraints = 2*(numSpatial) + 2*(numTemporal);
    
    if loadCSV
        Aineq = load('label_cache/Aineq.mat');
        Aineq = Aineq.Aineq;
        bineq = load('label_cache/bineq.mat');
        bineq = bineq.bineq;
        
        disp('Generating equality constraints');
        Aeq = load('label_cache/Aeq.mat');
        Aeq = Aeq.Aeq;
        beq = load('label_cache/beq.mat');
        beq = beq.beq; 
    else
        Aineq = zeros(numInequalityConstraints, numVariables);
        bineq = zeros(numInequalityConstraints, 1);
    end
    
    % helper functions for indexing into the arrays
    linearLabelIndex = @(i,j,t) ((j-1) + (i-1)*N + (t-1)*M*N + 1);
    linearSpatIndex = @(i,j,t,a,b) ((a-1) + (b-1)*spatialWinDim + ...
        (j-1)*spatialWinDim^2 + (i-1)*N*spatialWinDim^2 + ...
        (t-1)*M*N*spatialWinDim^2 + 1);
    linearTempIndex = @(i,j,t,a,b) ((a-1) + (b-1)*temporalWinDim + ...
        (j-1)*temporalWinDim^2 + (i-1)*N*temporalWinDim^2 + ...
        (t-1)*M*N*temporalWinDim^2 + 1);
    
    linearXIndex = @(i,j,t)     (startIndices(1) + linearLabelIndex(i,j,t) - 1);
    linearYIndex = @(i,j,t,a,b) (startIndices(2) + linearSpatIndex(i,j,t,a,b) - 1);
    linearZIndex = @(i,j,t,a,b) (startIndices(3) + linearTempIndex(i,j,t,a,b) - 1);
    
    if ~loadCSV
        % Constraint 1 - Y label similarity constraint
        disp('Generating inequality constraint 1');
        index = 1; % keeps track of the current constraint for ease of coding
        for t = 1:T
            for i = 1:M
                for j = 1:N
                    for b = 1:spatialWinDim
                        for a = 1:spatialWinDim
                            j_prime = safeIndex(j, a-spatialWinDim-1, N);
                            i_prime = safeIndex(i, b-spatialWinDim-1, M);

                            X_ijt_index = linearXIndex(i,j,t);
                            X_Y_index = linearXIndex(i_prime,j_prime,t);
                            Y_index = linearYIndex(i,j,t,a,b);

                            % add positive inequality constraint
                            Aineq(index+0, X_ijt_index) = 1;
                            Aineq(index+0, X_Y_index) = -1;
                            Aineq(index+0, Y_index) = -1;

                            % add negative inequality constraint
                            Aineq(index+1, X_ijt_index) = -1;
                            Aineq(index+1, X_Y_index) = 1;
                            Aineq(index+1, Y_index) = -1;

                            index = index+2;
                        end
                    end
                end
            end
        end

        % Constraint 2 - Z label similarity constraint
        disp('Generating inequality constraint 2');
        for t = 1:(T-1)
            for i = 1:M
                for j = 1:N
                    for b = 1:temporalWinDim
                        for a = 1:temporalWinDim
                            j_prime = safeIndex(j, a-temporalWinSize-1, N);
                            i_prime = safeIndex(i, b-temporalWinSize-1, M);

                            X_ijt_index = linearXIndex(i,j,t);
                            X_ijt_plus_index = linearXIndex(i_prime,j_prime,t+1);
                            Z_index = linearZIndex(i,j,t,a,b);

                            if i == 6 && j == 10 && a ==5 && b == 5
                               stop = 1; 
                            end
                            
                            % add positive inequality constraint
                            Aineq(index+0, X_ijt_index) = 1;
                            Aineq(index+0, X_ijt_plus_index) = -1;
                            Aineq(index+0, Z_index) = -1;

                            % add negative inequality constraint
                            Aineq(index+1, X_ijt_index) = -1;
                            Aineq(index+1, X_ijt_plus_index) = 1;
                            Aineq(index+1, Z_index) = -1;

                            index = index+2;
                        end
                    end
                end
            end
        end

        Aineq = sparse(Aineq);
        bineq = sparse(bineq);
        
        disp('Generating equality constraints');
        numEqualityConstraints = numPixelsPerFrame; % only constrain the first frame
        Aeq = zeros(numEqualityConstraints, numVariables);
        beq = zeros(numEqualityConstraints, 1);
        
        Aeq(1:numEqualityConstraints,1:numEqualityConstraints) = ...
            eye(numEqualityConstraints);
        k = 1;
        for i = 1:M
            for j = 1:N
                beq(k) = X{1}(i,j);
                k = k+1;
            end
        end
    end
    
    if saveCSV
        disp('Saving CSV files');
        save('label_cache/L.mat', 'L');
        save('label_cache/H.mat', 'H');
        save('label_cache/Aineq.mat', 'Aineq');
        save('label_cache/bineq.mat', 'bineq');
        save('label_cache/Aeq.mat', 'Aeq');
        save('label_cache/beq.mat', 'beq');
    end

    % Create the estimated flow vector
    X_init = zeros(numVariables, 1);
    startIndex = 1;
    endIndex  = startIndex + numPixelsPerFrame - 1;
    for t = 1:T
        X_init(startIndex:endIndex) = reshape(X{t}, [numPixelsPerFrame,1]);
    end
    
    disp('Solving linear program');
    [X_new, fval, exitflag, output] = cplexbilp(C, Aineq, bineq, Aeq, beq, X_init);

    if debug
        % query pixel index for debugging
        i = 10;
        j = 10;
        t = 1;
        a = 4;
        b = 4;
        
        % pull variable values for spcified index
        xi = linearXIndex(i,j,t);
        
        i_plus = i+b-temporalWinSize-1;
        j_plus = j+a-temporalWinSize-1;
        xi_plus = linearXIndex(i_plus,j_plus,t+1);
        
        neighborhood_x = zeros(spatialWinDim, spatialWinDim);
        neighborhood_y = zeros(spatialWinDim, spatialWinDim);
        for a = -spatialWinSize:spatialWinSize
            for b = -spatialWinSize:spatialWinSize
                i_neighbor = safeIndex(i, a, M);
                j_neighbor = safeIndex(j, b, N);
                x_index = linearXIndex(i_neighbor,j_neighbor,t+1);
                y_index = linearYIndex(i_neighbor,j_neighbor,t+1,a+spatialWinSize+1,b+spatialWinSize+1);
                
                neighborhood_x(a+spatialWinSize+1,b+spatialWinSize+1) = ...
                    X_new(x_index);
                neighborhood_y(a+spatialWinSize+1,b+spatialWinSize+1) = ...
                    X_new(y_index);
            end  
        end
        
        
        %yi = linearYIndex(i,j,t,a,b);
        zi = linearZIndex(i,j,t,a,b);
        x = X_new(xi);
        x_plus = X_new(xi_plus);
        %y = X_new(yi);
        z = X_new(zi);
        test = 1;
    end
          
    % recover flows cell matrix from output vector
    newX = cell(1,T);
    index = 1;
    for t = 1:T
        newX{t} = zeros(M, N); 
        for i = 1:M
            for j = 1:N
                % here we assume all t of U come before all t of V,
                % but above its U,V alternating!
                newX{t}(i, j) = X_new(startIndices(1)+index-1);
                index = index+1;
            end
        end
    end
end


