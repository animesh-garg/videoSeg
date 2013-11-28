function [ final_labels ] = propagate_labels( init_label, imSequence, W, params)
%PROPAGATE_FLOWS Solve for labels in next T time steps given all the
%weights.
%Input: 
%   init_label: the matrix containing the label for initial image (initial 
%       segmentation) 
%   imSequence: sequence of image of length T
%   T: number of time steps in future for prediction. (parameter)
%   W: weights 


lambda = params.lambda;
sigma = params.sigma;
spatial_nbd_size = params.spatial_nbd_size;
nbdSize = (2*spatial_nbd_size+1)^2;
%% Get cost acc. to appearance model for each frame t in T
[M,N,T]=size(imSequence); 
numLabels = T*M*N;

% gets costs
C  = getUnaryPotential(imSequence, init_label, threshEPS)

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
    for t = 1:T %Can do without this loop but makes it easier to understand
        cplex.addCols(lambda(1)*unaryCost(t,:),[],lb,ub,cType);        
    end
    
    
    numAuxVar_SLC = numLabels*(nbdSize^2); %this is BIG of the order of 2^25
    
    %adding auxiliary variables and constraints for spatial labelling coherence
    for t = 1:T              
        for i = 1: M
            for j = 1:N
                %Define Spatial Neighborhood: 
                % nbd_i = [i-spatial_nbd_size:i+spatial_nbd_size]
                for b = (-spatial_nbd_size):spatial_nbd_size
                    for a  = (-spatial_nbd_size):spatial_nbd_size                                                
                        nbd_i= min(max(i+a,1), M);
                        nbd_j = min(max(j+b,1), N);
                       
                        %add auxiliary variable for spatial labelling coherence
                        cplex.addCols(lambda(2),[], 0,1, 'C');
                        
                        coeffTemp1= sparse
                        cplex.addRows(-inf, coeffTemp1, 0);
                        cplex.addRows(-inf, coeffTemp2, 0);
                        
                    end
                end                
            end
        end
        
    end
    
    %adding auxiliary variables for temporal labelling coherence
    
    
        
    %Add constraints (try Parfor here for speedups)
    
    for t = 1:T             
        
        % add constraints for shrinkage/expansion
        
        %add constraints for temporal labelling coherence
        
    end

    cplex.solve()

catch m
    throw (m);

end

end

