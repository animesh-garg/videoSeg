function [ final_labels ] = propagate_labels( init_label, imSequence, W, params)
%PROPAGATE_FLOWS Solve for labels in next T time steps given all the
%weights.
%Input: 
%   init_label: the matrix containing the label for initial image (initial 
%       segmentation) 
%   imSequence: sequence of image of length T
%   T: number of time steps in future for prediction. (parameter)
%   W: weights 
%   params: parameters values for the current iteration

%% Get parameter values
lambda = params.lambda;
sigma = params.sigma;

spatial_nbd_size = params.spatial_nbd_size;
nbdSize_Spat = (2*spatial_nbd_size+1)^2;

temporal_nbd= params.h;
nbdSize_Temp = (2*temporal_nbd+1)^2;

if exists (params.threshEPS)
    threshEPS = params.threshEPS;
else
    threshEPS = 0.5 ; %Assign some useful Value ToDO
end

% If video sequence is multidim matrix
[M,N,T]=size(imSequence);

% if videosequence is a cell array
T = length(imSequence.I);
[M,N] = size(imSequence.I{1});


%% Get cost acc. to appearance model for each frame t in T 
% gets costs
C = getUnaryCosts(imSequence, init_label, threshEPS);

% rewrite Cost matrix C(T,M,N) into C(T, M*N)
unaryCost = reshape(C,T, M*N);

%% Get total number of variables in the forumulation
numLabels = T*M*N; %number of variables of the type X_ijt
numAuxVar_SLC = numLabels*(nbdSize_Spat); %this is BIG of the order of 2^25
numAuxVar_TLC = numLabels*(nbdSize_Temp); 

numVariables = numLabels+ numAuxVar_SLC + numAuxVar_TLC ;

%% Cplex Solver
try
    cplex = Cplex('propagate_labels');
    % Turn off CPLEX logging
    cplex.DisplayFunc = [];
    cplex.model.sense='minimize';
         
    %% Adding binary variables for X_ijt
    cType = char(ones(1, M*N)*'B');
    lb = zeros(1, M*N); ub = ones(1, M*N);   
    for t = 1:T %Can do without this loop but makes it easier to understand
        cplex.addCols(lambda(1)*unaryCost(t,:),[],lb,ub,cType);        
    end
    clear lb ub cType
    
    %% adding auxiliary variables and constraints for spatial labelling coherence
    varArray = [T,M,N];
    varAuxSpatArray = [T,M,N,2*spatial_nbd_size+1, 2*spatial_nbd_size+1];
    coeffTemp1= sparse(numVariables,1);
    coeffTemp2= sparse(numVariables,1);
    
    %Possible vectorization exists.
    % Make the idx_ijt as vectors and then then coeffTemp would be sparse
    % matrices instead of arrays.
    %ToDO: also if we define coeffTemp1 inside the loop, then we could
    %possibly use parfor or spmd
    for t = 1:T              
        for i = 1: M
            for j = 1:N
                %Define Spatial Neighborhood: 
                % nbd_i = [i-spatial_nbd_size:i+spatial_nbd_size]
                for a = (-spatial_nbd_size):spatial_nbd_size
                    for b  = (-spatial_nbd_size):spatial_nbd_size                                                
                        nbd_i= min(max(i+a,1), M);
                        nbd_j = min(max(j+b,1), N);
                       
                        %add auxiliary variable for spatial labelling coherence
                        cplex.addCols(lambda(2),[], 0,1, 'C');
                        
                        %add correspording pair of constraints                        
                        Idx_ijt = sub2ind(varArray, t,i,j); %index of x_ijt
                        IdX_ijt_y = sub2ind(varArray, t,nbd_i,nbd_j); %index of x_ijt_y
                        IdX_varAux = numLabels+ sub2ind(varAuxSpatArray, t,i,j,a,b);
                        
                        %Adding x_ijt - x_y <= p_ijty
                        coeffTemp1(Idx_ijt) = 1; 
                        coeffTemp1(Idx_ijt_y) = -1;
                        coeffTemp1(IdX_varAux) = -1;
                        cplex.addRows(-inf, coeffTemp1, 0);
                        
                        %reset coeffTemp1 for reuse in later iterations                        
                        coeffTemp1(Idx_ijt) = 0; 
                        coeffTemp1(Idx_ijt_y) = 0;
                        coeffTemp1(IdX_varAux) = 0;                        
                        
                        %Adding -x_ijt + x_y <= p_ijty
                        coeffTemp2(Idx_ijt) = -1; 
                        coeffTemp2(Idx_ijt_y) = 1;
                        coeffTemp2(IdX_varAux) = -1;
                        cplex.addRows(-inf, coeffTemp2, 0);                         
                        
                        %reset coeffTemp1 for reuse in later iterations                        
                        coeffTemp2(Idx_ijt) = 0; 
                        coeffTemp2(Idx_ijt_y) = 0;
                        coeffTemp2(IdX_varAux) = 0;
                    end
                end                
            end
        end        
    end
    
    clear coeffTemp1 coeffTemp2
    
    %% adding auxiliary variables and constraints for temporal labelling coherence        
    varAuxTempArray = [T,M,N,2*temporal_nbd+1, 2*temporal_nbd+1];
    coeffTemp1= sparse(numVariables,1);
    coeffTemp2= sparse(numVariables,1);
    
    for t = 1:T-1              
        for i = 1: M
            for j = 1:N
                %Define Temporal Neighborhood: 
                % nbd_i = [i-temporal_nbd:i+temporal_nbd]
                for a = (-temporal_nbd):temporal_nbd
                    for b  = (-temporal_nbd):temporal_nbd                                                
                        nbd_i= min(max(i+a,1), M);
                        nbd_j = min(max(j+b,1), N);
                       
                        %add auxiliary variable for temporal labelling coherence
                        cplex.addCols(lambda(3)*W{t}(i,j,a,b),[], 0,1, 'C');                        
                        
                        %add correspording pair of constraints                        
                        Idx_ijt = sub2ind(varArray, t,i,j); %index of x_ijt
                        %index of x_ijtPlus_ab
                        IdX_ijtPlus_ab = sub2ind(varArray, (t+1),nbd_i,nbd_j); 
                        IdX_varAux = numLabels+ numAuxVar_SLC +...
                            sub2ind(varAuxSpatArray, t,i,j,a,b);
                        
                        %Adding x_ijt - x_ijtPlus_ab <= q_ijtab
                        coeffTemp1(Idx_ijt) = 1; 
                        coeffTemp1(IdX_ijtPlus_ab) = -1;
                        coeffTemp1(IdX_varAux) = -1;
                        cplex.addRows(-inf, coeffTemp1, 0);
                        
                        %reset coeffTemp1 for reuse in later iterations                        
                        coeffTemp1(Idx_ijt) = 0; 
                        coeffTemp1(IdX_ijtPlus_ab) = 0;
                        coeffTemp1(IdX_varAux) = 0;                        
                        
                        %Adding -x_ijt + x_ijtPlus_ab <= q_ijtab
                        coeffTemp2(Idx_ijt) = -1; 
                        coeffTemp2(IdX_ijtPlus_ab) = 1;
                        coeffTemp2(IdX_varAux) = -1;
                        cplex.addRows(-inf, coeffTemp2, 0);                         
                        
                        %reset coeffTemp2 for reuse in later iterations
                        coeffTemp2(Idx_ijt) = 0; 
                        coeffTemp2(IdX_ijtPlus_ab) = 0;
                        coeffTemp2(IdX_varAux) = 0;
                    end
                end                
            end
        end        
    end
    
    clear coeffTemp1 coeffTemp2
    
    %% add constraints for shrinkage/expansion
    coeffTemp1= sparse(numVariables,1);
    coeffTemp2= sparse(numVariables,1);
    [I,J] = ind2sub([M,N],1:(M*N)); %generate indices for all pixels in one image
    for t = 1:T-1             
        %vector of indices for label variables in frame t 
        Idx_ijt=sub2ind(varArray, t*ones(M*N, 1),I',J'); 
        %vector of indices for labels in frame t+1
        Idx_ijtPlus=sub2ind(varArray, (t+1)*ones(M*N, 1),I',J');
        
        % add constraints FGinFrameT - FGinFrameTplus <= sigma*FGinFrameT 
        coeffTemp1 (Idx_ijt) = 1-sigma;
        coeffTemp1 (Idx_ijtPlus) = -1;
        cplex.addRows(-inf, coeffTemp1, 0);                         
        %reset coeffTemp1 for reuse in later iterations
        coeffTemp1 (Idx_ijt) = 0;
        coeffTemp1 (Idx_ijtPlus) = 0;
        
        % add constraints FGinFrameTplus - FGinFrameT <= sigma*FGinFrameT 
        coeffTemp2 (Idx_ijt) = (1+sigma);
        coeffTemp2 (Idx_ijtPlus) = 1;
        cplex.addRows(-inf, coeffTemp2, 0);
        %reset coeffTemp2 for reuse in later iterations
        coeffTemp2 (Idx_ijt) = 0;
        coeffTemp2 (Idx_ijtPlus) = 0;
        
    end
    clear coeffTemp1 coeffTemp2
    
    %% Add constraint x_1ij = init_label
    %There are other better ways (do this while defining objective but it 
    % hampers understanding    
    coeffTemp= sparse(numVariables,1);
    
    %Using I and J from previous set of constraints
    % I,J are indices for all pixels in one image as a list of tuples
    %vector of indices for label variables in frame 1 
    Idx_ijt=sub2ind(varArray, ones(M*N, 1),I',J');
    
    %add (M*N) equality constraints
    for k = 1: (M*N)
        coeffTemp(Idx_ijt(k)) =1;
        % if init_image pixel is FG, then add a constraint making the
        % corresponding label variable to be 1
        if (init_image(I(k), J(k))==1)            
            cplex.addRows(1, coeffTemp,1);
        % if init_image pixel is BG, then add a constraint making the
        % corresponding label variable to be 0
        elseif (init_image(I(k), J(k))==0)                  
            cplex.addRows(0, coeffTemp,0);
        else
            fprintf('Initial Mask needs to be Binary')
            return;
        end
        coeffTemp(Idx_ijt(k)) =0;
    end                
    
    clear coeffTemp
        
%% Call Cplex Solver
    cplex.solve()
    
    
    % Write the solution
    fprintf('\nSolution status = %s\n',cplex.Solution.statusstring);
    fprintf('\nSolution Time = %d\n', cplex.Solution.time);
    fprintf('Solution value = %f\n',cplex.Solution.objval);
    
    opt_sol = cplex.Solution.x;        

catch m
    throw (m);
    disp(m.message);

%% read results
% only the first T*M*N variable are the labels
final_labels = cell (1,T);

%final_labels{1} = reshape(opt_sol(1:M*N)', M, N); 
for t = 0: T-1
   final_labels{t+1} = reshape(opt_sol ((t*M*N +1) : (t+1)*M*N)', M,N);  
end
    
    
end

end
