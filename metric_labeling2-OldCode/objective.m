%Assignment of objective function  coefficients
function obj_coeff = objective(D, C, width, height,nbd_size )
    fprintf('starting objective function coefficient assignment\n');
    num_pixels= width*height;
    num_nbd=nbd_size;
    obj_coeff = sparse(2*num_pixels+(2*(num_pixels*num_nbd)),1);
    percent_complete = 0; 
    for i = 1:num_pixels 
        %defining coefficients of objective function
        %for x_i_0
        obj_coeff(i,1)=C(i,1)-C(i,2)+ 1.0000e-012;
        %for x_i_1
        obj_coeff(num_pixels+i,1)=0+ 1.0000e-012;
        
        %If corner pixel nbd_i=4
        if (i == 1||i==width || i== ((height-1)*width)+1|| i== num_pixels)
            nbd_i=4;
            if(i==1)
                nbd= [0 1 width width+1];
            elseif(i==width)
                nbd = [-1 0 width-1 width];
            elseif(i == ((height-1)*width)+1)
                nbd= [-width -(width-1) 0 1];
            else
                nbd= [-(width+1) -width -1 0];
            end
            
            for j = 1: nbd_i
                %for y01_i_j           
                obj_coeff(2*num_pixels+((i-1)*num_nbd)+j,1)=D(i,i+nbd(j))+ 1.0000e-012;
                %for y10_i_j
                obj_coeff((2*num_pixels)+(num_pixels*num_nbd)+(((i-1)*num_nbd)+j),1)=D(i,i+nbd(j))+ 1.0000e-012;
            end
            
          
            
        %if top or bottom edge pixel , nbd_i= 6
        elseif (i >1 && i<width)|| (i>((height-1)*width)+1 && i< num_pixels)
            nbd_i=6;
            if  (i >1 && i<width)
                nbd= [-1 0 1 width-1 width width+1];
            elseif (i>((height-1)*width)+1 && i< num_pixels)
                nbd = [-(width+1) -width -(width-1) -1 0 1];
            end
              
            for j = 1: nbd_i
                %for y01_i_j           
                obj_coeff(2*num_pixels+(((i-1)*num_nbd)+j),1)=D(i,i+nbd(j))+ 1.0000e-012;
                %for y10_i_j
                obj_coeff((2*num_pixels)+(num_pixels*num_nbd)+(((i-1)*num_nbd)+j),1)=D(i,i+nbd(j))+ 1.0000e-012;
            end
            
        %if left or right edge pixel , nbd_i= 6    
        elseif ( (mod(i,width)==1 ) &&  (i~=1 || i~=(height-1)*width+1 ) ) || ( (mod (i,width)==0) && (i~= width || i ~= num_pixels) ) 
                nbd_i = 6;
                if ( (mod(i,width)==1 ) &&  (i~=1 || i~=(height-1)*width+1 ) )
                    nbd= [-width -(width-1) 0 1 width width+1];
                elseif ( (mod (i,width)==0) && (i~= width || i ~= num_pixels) )
                    nbd = [-(width+1) -width -1 0 width-1 width];
                end
                for j = 1: nbd_i
                %for y01_i_j           
                obj_coeff(2*num_pixels+((i-1)*num_nbd)+j,1)=D(i,i+nbd(j))+ 1.0000e-012;
                %for y10_i_j
                obj_coeff((2*num_pixels)+(num_pixels*num_nbd)+(((i-1)*num_nbd)+j),1)=D(i,i+nbd(j))+ 1.0000e-012;
                end
            
            
        %Else every pixel in the middle nbd_i=9
        else
            nbd_i=9;
            nbd= [-(width+1) -width -(width-1) -1 0 1 width-1 width width+1];
            for j = 1: nbd_i
                %for y01_i_j           
                obj_coeff(2*num_pixels+((i-1)*num_nbd)+j,1)=D(i,i+nbd(j))+ 1.0000e-012;
                %for y10_i_j
                obj_coeff((2*num_pixels)+(num_pixels*num_nbd)+(((i-1)*num_nbd)+j),1)=D(i,i+nbd(j))+ 1.0000e-012;
            end
        end
        
        %output the rate of completion
        if (round(100*i/num_pixels)==percent_complete+20)
            percent_complete = percent_complete+20;
            fprintf('Processing %02i%% \n',percent_complete)
        end
    end