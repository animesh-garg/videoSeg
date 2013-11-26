<<<<<<< HEAD
function [videoStruct] = read_data(videoName)
%READ_DATA Returns a struct for the given video name
%   Input : videoName (eg - 'Cars1')
%   Output : VideoStruct with fields I (frame sequence), X1(initial segmentation)


%% Initializing Paths
startup;
videoPath = sprintf('%s%s', datasetPath, videoName);

%% Getting names of image files
imageNames = getFileNamesFromDirectory(videoPath,'mode','path','types',{'.jpg','.jpeg','.png'}); %Some script i'd written earlier. Included in utils
T = length(imageNames);

%% Reading the frames of the video
I = {};
for t=1:T
    I{t} = imresize(imread(imageNames{t}),d); %read and downsample image
end

%% Obtaining Initial segmentation for some object in the video
segNames = getFileNamesFromDirectory(videoPath,'mode','path','types',{'01.pgm'});
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
=======
function [ output_args ] = read_data( input_args )
%READ_DATA Summary of this function goes here
%   Detailed explanation goes here


end

>>>>>>> 78812a5e33655e88215d3ad9256cba61958cf9a4
