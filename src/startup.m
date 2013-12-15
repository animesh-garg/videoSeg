addpath(genpath('utils'));
addpath(genpath('AppearanceCostFunctions'));
addpath(genpath('solvers'));
datasetPath = '../data/SegTrackv2/JPEGImages/';
groundTruthPath = '../data/SegTrackv2/GroundTruth/';

b = 2; %Neighbourhood Size for next frame in which current pixel would map to
d = 1; %Downsample rate