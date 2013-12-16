function [videoStruct] = read_data(videoName, T, startFrame)
%READ_DATA Returns a struct for the given video name
%   Input : videoName (eg - 'Cars1')
%   Output : VideoStruct with fields I (frame sequence), X1(initial segmentation)


%% Initializing Paths
startup;
videoPath = sprintf('%s%s', datasetPath, videoName);
gtPath = sprintf('%s%s', groundTruthPath, videoName);

%% Getting names of image files
imageNames = getFileNamesFromDirectory(videoPath,'mode','path','types',{'.jpg','.jpeg','.png','.bmp'}); %Some script i'd written earlier. Included in utils
if nargin < 2
    T = length(imageNames);
end
if nargin < 3
    startFrame = 1;
end

%% Reading the frames of the video
I = {};
for t=startFrame:(startFrame+T)
    I{t-startFrame+1} = imresize(imread(imageNames{t}),d); %read and downsample image
end

%% Obtaining Initial segmentation for some object in the video
segNames = getFileNamesFromDirectory(gtPath,'mode','path','types',{'15.png'});
if (~isempty(segNames))
    X1 = imresize(imread(segNames{1}),d);
    %initialize for one of the available segmented objects
    X1(X1 < 255) = 0;
    X1(X1==255) = 1;
else
    disp('Error: No segments available')
    X1 = [];
end

%% Returning the videoStruct
videoStruct = struct();
videoStruct.I = I;
videoStruct.X1 = X1;

end
