addpath(genpath('utils'));
addpath(genpath('AppearanceCostFunctions'));
addpath(genpath('solvers'));
datasetPath = '../data/moseg_dataset/';

b = 2; %Neighbourhood Size for next frame in which current pixel would map to
d = 1; %Downsample rate