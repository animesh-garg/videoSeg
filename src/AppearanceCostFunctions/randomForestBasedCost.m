function [cost] = randomForestBasedCost(vid)
%DECISIONTREEBASEDCOST Gives foreground cost for each pixel
%   Input : VideoStruct
%   Output : cell of length T. Each cell element is a height X widhth
%   matrix with value = foreground cost for that pixel (between [0 1]). The
%   background cost is 1 - foreground Cost
%   Even though this is slightly slower, I would suggest use of this over average value based model as it is
%   more reliable as it handles variations within foreground and background well

%% Initializing Variables

I = vid.I;
X1 = vid.X1;
w = 2; %patch Width = 2*w+1 X 2*w+1
T = length(I);
[ydim,xdim,C] = size(I{1});
frame1 = zeros(ydim+2*w,xdim+2*w,C);
frame1(w+1:ydim+w,w+1:xdim+w,:) = double(I{1});

%% Training Data for classifier
pTrain={'maxDepth',25,'M',20,'minChild',5};
labels = 1+uint8(X1(:));
data = zeros(ydim*xdim,C*(2*w+1)^2);
index = 0;
for x = w+1:xdim+w
    for y = w+1:ydim+w
        index = index+1;
        patch = frame1(y-w:y+w,x-w:x+w,:);
        data(index,:) = patch(:)';
    end
end

%% Training Classifier
forest = forestTrain(data, labels, pTrain{:});
disp('Trained Classfier')

%% Computing Costs
cost = cell(T,1);
cost{1} = 2*(ones(ydim,xdim)-double(X1))-1;

for t=2:T
    fprintf ('Processing frame number %d of %d \r', t,T);
    frame = zeros(ydim+2*w,xdim+2*w,C);
    frame(w+1:ydim+w,w+1:xdim+w,:) = double(I{t});
    data = zeros(ydim*xdim,C*(2*w+1)^2);
    
    index = 0;
    for x = w+1:xdim+w
        for y = w+1:ydim+w
            index = index+1;
            patch = frame(y-w:y+w,x-w:x+w,:);
            data(index,:) = patch(:)';
        end
    end
    
    [~,probs] = forestApply(single(data),forest);
    cost{t} = reshape(probs(:,1),[ydim,xdim]);
    cost{t} = 2*cost{t} - 1;
    
end


end

