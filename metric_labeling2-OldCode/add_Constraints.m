function [cplex, t_2]= add_Constraints(num_pixels,nbd_size, cplex )
    %variable counter for completion
    percent_complete=0;

    time_constraint_assignment = tic;
    
    %Adding constraints one by one
    %y01_i_j-x_i_0<=0           
    for i = 1:num_pixels
        for j = 1: nbd_size
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A(i,1) = -1;
            A((2*num_pixels)+((i-1)*nbd_size+j),1)=1;
            cplex.addRows(-1,A(:,1)' ,0); 
            %clear A;
        end
        %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint y01_i_j-x_i_0<=0.. %02i%% \n',percent_complete)
        end
    end
    
    percent_complete=0;
    %y01_i_j-x_j_1<=0
    for i = 1:num_pixels
        for j = 1: nbd_size
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A(num_pixels+j,1)=-1;
            A((2*num_pixels)+(((i-1)*nbd_size)+j),1)=1;
            cplex.addRows(-1,A(:,1)' ,0); 
        end
         %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint y01_i_j-x_j_1<=0.. %02i%% \n',percent_complete)
        end
    end
    
    percent_complete=0;
    %y10_i_j-x_i_1<=0
    for i = 1:num_pixels
        for j = 1: nbd_size
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A(num_pixels+i,1)=-1;
            A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size)+j),1)=1;
            cplex.addRows(-1,A(:,1)' ,0); 
        end
         %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint y10_i_j-x_i_1<=0.. %02i%% \n',percent_complete)
        end
    end
    
    percent_complete=0;
    %y10_i_j-x_j_0<=0
    for i = 1:num_pixels
        for j = 1: nbd_size
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A(j,1)=-1;
            A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size+j)),1)=1;
            cplex.addRows(-1,A(:,1)' ,0); 
        end
        %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint y10_i_j-x_j_0<=0.. %02i%% \n',percent_complete)
        end 
    end
    
    percent_complete=0;
    %x_i_0+x_j_1-y01_i_j<=1
    for i = 1:num_pixels
        for j = 1: nbd_size
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A(i,1)=1;
            A(num_pixels+j,1)=1;
            A((2*num_pixels)+((i-1)*nbd_size+j),1)=-1;
            cplex.addRows(0,A(:,1)' ,1); 
        end
         %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint x_i_0+x_j_1-y01_i_j<=1.. %02i%% \n',percent_complete)
        end
    end
    
    percent_complete=0;
    %x_i_1+x_j_0-y10_i_j<=1
    for i = 1:num_pixels
        for j = 1: nbd_size
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A(num_pixels+i,1)=1;
            A(j,1)=1;
            A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size+j)),1)=-1;
            cplex.addRows(0,A(:,1)' ,1); 
        end
        %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint x_i_1+x_j_0-y10_i_j<=1.. %02i%% \n',percent_complete)
        end 
    end
    
    percent_complete=0;
    %0<=y01_i_j+y10_i_j<=1
    for i = 1:num_pixels
        for j = 1: nbd_size
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A((2*num_pixels)+(((i-1)*nbd_size+j)),1)=1;
            A((2*num_pixels)+(num_pixels*nbd_size)+(((i-1)*nbd_size+j)),1)=1;
            cplex.addRows(0,A(:,1)' ,1); 
        end
        %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint 0<=y01_i_j+y10_i_j<=1.. %02i%% \n',percent_complete)
        end 
    end
    
    percent_complete=0;
    %1<=x_i_0+x_i_1<=1
    for i = 1:num_pixels
            A = sparse(2*num_pixels+(2*(num_pixels*nbd_size)),1);
            A(i,1)=1;
            A(num_pixels+i,1)=1;
            cplex.addRows(1,A(:,1)' ,1); 
            %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing constraint 1<=x_i_0+x_i_1<=1.. %02i%% \n',percent_complete)
        end 
    end
     
    t_2= toc(time_constraint_assignment);