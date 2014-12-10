function [Cost] = generate_imCost(init_image)
    [height, width] = size(init_image);
    w_data=0.5;
    w_smoothness=0.5;
    padded_image = padarray (init_image, [1 1]); 
  
   for i= 2:width
        for j =2:height
            count_zeros_in_nbd = 0;
            count_ones_in_nbd = 0;
            if (init_image(j,i)==0)
                for l = i-1:i+1
                    for k = j-1:j+1
                        if (padded_image(k,l) == 0)
                            count_zeros_in_nbd = count_zeros_in_nbd+ 1;
                        end
                    end
                end
                %update cost matrix
                % the cost for each term is normalized
                %notice that ((i-1)*j+j)= pixel_idx
                nbd_Cost(((i-1)*j+j),1)=(count_zeros_in_nbd-1)/9;
            else
                 for l = i-1:i+1
                    for k = j-1:j+1
                        if (padded_image(k,l) == 1)
                            count_ones_in_nbd = count_ones_in_nbd+1;
                        end
                    end
                end
                %update cost matrix
                % the cost for each term is normalized
                %notice that ((i-1)*j+j)= pixel_idx
                nbd_Cost(((i-1)*j+j),2)=(count_ones_in_nbd-1)/9;
            end
        end
   end
   
   data_cost = zeros(width*height,2);
   for i= 1:width
        for j =1:height
            if(init_image(j,i)==1)
                data_cost((i-1)*j+j,1) =1;
                data_cost((i-1)*j+j,2) =0;
            else
                data_cost((i-1)*j+j,1) =0;
                data_cost((i-1)*j+j,2) =1;
            end
         end
   end
   
   %the smoothness term is subtracted because nbd_Cost is a variable which
   %counts number of neigghtbors of same attribute.
   Cost = w_data*data_cost - w_smoothness*nbd_Cost;
end

