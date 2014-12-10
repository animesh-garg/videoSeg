function [D] = generate_dist(init_image)
    [height,width] = size(init_image);
        
    %assigning every pixel a pixel_idx = (i-1)*j+j where k and l are loop
    %counters in the current iteration
    %Debugging: This allocation is Giving memory overflow error. hence we deploy write
    %to .mat file option
    %D = zeros (width*height,width*height);
    %save ('dist.mat','init_image')
    D = sparse(width*height, width*height);
    %top left corner pixel
    D (1, 2) = 1;
    D (1, width+1) = 1;
    D (1, width +2) = 2;
    %top right corner pixel
    D (width, width-1) = 1;
    D (width, 2*width) = 1;
    D (width, 2*width -1) = 2;
    %bottom left corner
    D ((height-1)*width+1, (height-1)*width+2) = 1;
    D ((height-1)*width+1, (height-2)*width+1) = 1;
    D ((height-1)*width+1, (height-2)*width+2) = 2;
    %bottom right corner
    D (width*height, width*height-1) = 1;
    D (width*height, width*(height-1)) = 1;
    D (width*height, width*(height-1)-1) = 2;
    
     
    for i = 2:width-1
        %top edge
        D(i,i-1)=1;
        D(i,i+1)=1;
        D(i,width-1)=2;
        D(i, width)=1;
        D(i,width+1)=2;
        %bottom edge
        last_row = (height-1)*width;
        D(last_row+i,last_row+i-1)=1;
        D(last_row+i,last_row+i+1)=1;
        D(last_row+i,last_row-width+1)=2;
        D(last_row+i,last_row-width+2)=1;
        D(last_row+i,last_row-width+3)=2;
    end
    fprintf('\ntop and bottom edges completed...\n');
    for i = 2: height-1
        %left edge
        D(width*(i-1)+1,width*(i-2)+1)=1;
        D(width*(i-1)+1,width*(i-2)+2)=2;
        D(width*(i-1)+1,width*(i-1)+2)=1;
        D(width*(i-1)+1,(width*i)+1)=1;
        D(width*(i-1)+1,(width*i)+2)=2;
        %right edge
        D(width*i, width)=1;
        D(width*i, width-1)=2;
        D(width*i, (width*i) - 1)=1;
        D(width*i, width*(i+1))=1;
        D(width*i, width*(i+1)-1)=1;
    end
    fprintf('left and right edges completed...\n');
    for i = 2 : height-1
        for j = 2: width -1
            D(i*width+j, (i-1)*width+(j-1))=2; 
            D(i*width+j, (i-1)*width+j)=1;
            D(i*width+j, (i-1)*width+(j+1))=2;
            D(i*width+j, i*width+(j-1))=1;
            D(i*width+j, i*width+(j+1))=1;
            D(i*width+j, (i+1)*width+(j-1))=2;
            D(i*width+j, (i+1)*width+j)=1;
            D(i*width+j, (i+1)*width+(j+1))=2;
        end
       % save ('dist.mat', 'D','-append')
    end
    fprintf('all of Distances completed\n');
end
