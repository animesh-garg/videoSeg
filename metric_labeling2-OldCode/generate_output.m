function generate_output()
    addpath ('C:\ILOG\CPLEX121\matlab\x64_win64')
try
    time_count= tic;
    path = './test_bw_05.tif';
    init_image = read_image(path);
    nbd_size = 9;
    D =generate_dist(init_image);
    C= generate_imCost(init_image);
    
    [height, width] = size(init_image);
    num_pixels= width*height;
   
    
    cplex= Cplex ('generate_output');
    cplex.readModel('metric_labeling.lp');
    cplex = add_objective(D, C, width, height, nbd_size, cplex);
    cplex.writeModel('metric_labeling2.lp');
    cplex.solve();
    
    
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
    
    t1 =toc(time_count);
    toc(time_count)
catch m
     disp(m.message);
end