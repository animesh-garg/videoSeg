function ML ()
    clear all;
    addpath ('C:\ILOG\CPLEX121\matlab\x64_win64');
try
    %counter for time total
    time_count= tic;
    path = './test_bw_05.tif';
    init_image = read_image(path);
    D =generate_dist(init_image);
    C= generate_imCost(init_image);
    %load ('dist.mat','D')
    [height, width] = size(init_image);
    num_pixels= width*height;
    %size of neighborhood for computation (include the pixel itself)
    %dont change will affect other calculations
    nbd_size=9;
    %specify lower nd upper bounds
    lb    = sparse(2*num_pixels+2*(num_pixels*nbd_size),1);
    ub    = (1:(2*num_pixels+2*(num_pixels*nbd_size)))';
    
    % Variables types used for the Binary IP formulation
    x_i= char(ones([1 2*num_pixels])*'B');
    y_i_j = char(ones([1 2*num_pixels*nbd_size])*'B');
    
    % Specify variable types when integrality needed
    temp = strcat(x_i,y_i_j);
    clear x_i;
    clear y_i_j;
    ctype = char (temp);
    clear temp;
    
    %populate objective function coefficiant Column Matrix
    obj_coeff = objective (D, C, width, height,nbd_size);
     
    % Initialize the CPLEX object
    cplex = Cplex('ML');
    
    % Turn off CPLEX logging
    cplex.DisplayFunc = [];
    
    % Use addRows to populate the model
    % This is for adding just the objective function
    cplex.addCols(obj_coeff, [], lb, ub,ctype);
    cplex.Model.sense = 'minimize';
          
    %add constraints
    [cplex, t_2]= add_Constraints(num_pixels,nbd_size, cplex);

    %write model
    cplex.writeModel('metric_labeling_before_solving.lp');
    
    % Solve model
    cplex.solve();
    
    %write model
    cplex.writeModel('metric_labeling.lp');
    
    % Write the solution
    fprintf('\nSolution status = %s\n',cplex.Solution.statusstring);
    fprintf('\nSolution Time = %d\n', cplex.Solution.time);
    fprintf('Solution value = %f\n',cplex.Solution.objval);
    
    opt_sol = cplex.Solution.x;
    %make an image from solution
     temp = zeros (height, width);
    for i = 1: num_pixels 
        if ((mod(i,width)~=0) && (opt_sol(i,1)==1))
            temp(1+floor(i/width), mod(i,width))=0;
        elseif (mod(i,width)~=0 && opt_sol(num_pixels+i,1)==1)
            temp(1+floor(i/width), mod(i,width))=1;
        elseif ((mod(i,width) == 0) && (opt_sol(i,1) ==1))
            temp (1+floor(i/width), width)=0;
        elseif ((mod(i,width) == 0) && (opt_sol(num_pixels+i,1) ==1))
            temp (1+floor(i/width), width)=1;
        else   
            fprintf('Error in Solution. Neither label assigned to a pixelid (%d)\n',i);
        end
    end
    
    s2=sprintf('./test_bw_05_result.tif');
    imwrite(temp,s2);
    
    t_1 =toc(time_count);
    disp('t_2=',t_2);
    disp ('t1=',t_1);
catch m
     disp(m.message);
end




%     %variable counter for completion
%     percent_complete=0;
%     
%     time_constraint_assignment = tic;
%     %Adding constraints one by one
%     %y01_i_j-x_i_0<=0           
%     for i = 1:num_pixels
%         for j = 1: nbd_size
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A(i,1) = -1;
%             A((2*num_pixels)+((i-1)*nbd_size+j),1)=1;
%             cplex.addRows(-1,A(:,1)' ,0); 
%             %clear A;
%         end
%         %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint y01_i_j-x_i_0<=0.. %02i%% \n',percent_complete)
%         end
%     end
%     
%     percent_complete=0;
%     %y01_i_j-x_j_1<=0
%     for i = 1:num_pixels
%         for j = 1: nbd_size
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A(num_pixels+j,1)=-1;
%             A((2*num_pixels)+(((i-1)*nbd_size)+j),1)=1;
%             cplex.addRows(-1,A(:,1)' ,0); 
%         end
%          %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint y01_i_j-x_j_1<=0.. %02i%% \n',percent_complete)
%         end
%     end
%     
%     percent_complete=0;
%     %y10_i_j-x_i_1<=0
%     for i = 1:num_pixels
%         for j = 1: nbd_size
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A(num_pixels+i,1)=-1;
%             A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size)+j),1)=1;
%             cplex.addRows(-1,A(:,1)' ,0); 
%         end
%          %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint y10_i_j-x_i_1<=0.. %02i%% \n',percent_complete)
%         end
%     end
%     
%     percent_complete=0;
%     %y10_i_j-x_j_0<=0
%     for i = 1:num_pixels
%         for j = 1: nbd_size
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A(j,1)=-1;
%             A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size+j)),1)=1;
%             cplex.addRows(-1,A(:,1)' ,0); 
%         end
%         %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint y10_i_j-x_j_0<=0.. %02i%% \n',percent_complete)
%         end 
%     end
%     
%     percent_complete=0;
%     %x_i_0+x_j_1-y01_i_j<=1
%     for i = 1:num_pixels
%         for j = 1: nbd_size
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A(i,1)=1;
%             A(num_pixels+j,1)=1;
%             A((2*num_pixels)+((i-1)*nbd_size+j),1)=-1;
%             cplex.addRows(0,A(:,1)' ,1); 
%         end
%          %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint x_i_0+x_j_1-y01_i_j<=1.. %02i%% \n',percent_complete)
%         end
%     end
%     
%     percent_complete=0;
%     %x_i_1+x_j_0-y10_i_j<=1
%     for i = 1:num_pixels
%         for j = 1: nbd_size
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A(num_pixels+i,1)=1;
%             A(j,1)=1;
%             A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size+j)),1)=-1;
%             cplex.addRows(0,A(:,1)' ,1); 
%         end
%         %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint x_i_1+x_j_0-y10_i_j<=1.. %02i%% \n',percent_complete)
%         end 
%     end
%     
%     percent_complete=0;
%     %0<=y01_i_j+y10_i_j<=1
%     for i = 1:num_pixels
%         for j = 1: nbd_size
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A((2*num_pixels)+(((i-1)*nbd_size+j)),1)=1;
%             A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size+j)),1)=1;
%             cplex.addRows(0,A(:,1)' ,1); 
%         end
%         %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint 0<=y01_i_j+y10_i_j<=1.. %02i%% \n',percent_complete)
%         end 
%     end
%     
%     percent_complete=0;
%     %1<=x_i_0+x_i_1<=1
%     for i = 1:num_pixels
%             A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
%             A(i,1)=1;
%             A(num_pixels+i,1)=1;
%             cplex.addRows(1,A(:,1)' ,1); 
%             %output the rate of completion
%         if (round(100*i/num_pixels)==percent_complete+20)
%             percent_complete = percent_complete+20;
%             fprintf('Processing constraint 1<=x_i_0+x_i_1<=1.. %02i%% \n',percent_complete)
%         end 
%     end
%      
%     t_2= toc(time_constraint_assignment);