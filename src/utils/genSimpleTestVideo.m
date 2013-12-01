function [videoStruct] = genSimpleTestVideo()
    T = 3;
    m = 7;
    n = 7;
    % generate simple test video in grayscale
    
    videoStruct = struct;
    videoStruct.T = T;
    
    I = cell(1,T);
    I{1} = zeros(m,n);
    I{1}(2:3,2:3) = 255*ones(2,2);
    
    I{2} = zeros(m,n);
    I{2}(3:4,3:4) = 255*ones(2,2);
    
    I{3} = zeros(m,n);
    I{3}(4:5,4:5) = 255*ones(2,2);
    
    X1 = zeros(m,n);
    X1(2:3,2:3) = ones(2,2);
    
    videoStruct.I = I;
    videoStruct.X1 = X1;
end

