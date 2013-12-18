addpath(genpath('utils'));
addpath(genpath('AppearanceCostFunctions'));
addpath(genpath('solvers'));

addpath (genpath('classic_NL_code')) % to use the new optical flow solver


addpath (genpath('label_cache'))
addpath (genpath('flow_cache'))
addpath (genpath('weight_cache'))

datasetPath = '../data/moseg_dataset/';
groundTruthPath = '../data/moseg_dataset/';



b = 2; %Temporal Neighbourhood Size for next frame in which current pixel would map to
d = 1; %Downsample rate