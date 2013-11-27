function [ final_labels ] = propagate_labels( init_label, imSequence, W)
%PROPAGATE_FLOWS Solve for labels in next T time steps given all the
%weights.
%Input: 
%   init_label: the matrix containing the label for initial image (initial 
%       segmentation) 
%   imSequence: sequence of image of length T
%   T: number of time steps in future for prediction. (parameter)
%   W: weights 
    
%% Get cost acc. to appearance model for each frame t in T
[M,N,T]=size(imSequence); 
numLabels = T*M*N;

% gets costs
C  = getUnaryPotential(imSequence, init_label)

% rewrite Cost matrix C(T,M,N) into C(T, M*N)
unaryCost = reshape(C,T, M*N);

%% 



%% Cplex Solver
try
    cplex = Cplex('propagate_labels');
    % Turn off CPLEX logging
    cplex.DisplayFunc = [];
    cplex.model.sense='minimize';
    
    %add variables 
    %Adding binary variables for X_ijt
    cType = char(ones(1, M*N)*'B');
    lb = zeros(1, M*N); ub = ones(1, M*N);   
    for t = 1:T
        cplex.addCols(unaryCost(t,:),[],lb,ub,cType)
        
    end
    
    %adding auxiliary variables Y for the neighborhood constraints
    
    
    
    %add constraints (try Parfor here for speedups)
    
    for t = 1:T             
        %add constraints for spatial neighbourhood similarity
        for i = 1: M
            for j = 1:N
                %Define Neighborhood
                
            end
        end
        
        % add constraints for shrinkage/expansion
        
    end




    cplex.solve()
catch m
    throw (m);

end

end

