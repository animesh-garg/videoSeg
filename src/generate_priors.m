% Script to generate priors for our video segementatino algorithm
%   Uses the moseg_dataset and 3rd party NCuts code
%   Author: Jeff 

clear; clc; close all;

%% Get all images in a video sequence

b = 2;
nbSegments = 10;
dir_contents = dir('moseg_dataset');
video_folders = {};

n = 1;
for i = 1:size(dir_contents,1)
    if dir_contents(i).isdir && ~strcmp(dir_contents(i).name, '.') && ...
            ~strcmp(dir_contents(i).name, '..')
        video_folders{n} = dir_contents(i).name;
        n = n+1;
    end
end
num_videos = size(video_folders,2);

% now perform preprocessing for each folder
for i = 1:1 %num_videos
    video_folder = video_folders{i};
    
    % meat of the algo
    folder_files = dir(sprintf('moseg_dataset/%s', video_folder));
    video_frames = {};
   
    % read in all of the frames
    k = 1;
    for j = 1:size(folder_files,1)
        name = folder_files(j).name;
        if length(name) > 4 && ...
                strcmp(name((size(name,2)-3):size(name,2)), '.jpg')
        
            name = sprintf('moseg_dataset/%s/%s',video_folder, name);
            fprintf('%s\n',name);
            video_frames{k} = imread(name);
            k = k+1;
        elseif length(name) > 6 && ...
                strcmp(name((size(name,2)-5):size(name,2)), '01.pgm')
            name = sprintf('moseg_dataset/%s/%s',video_folder, name);
            initial_segment = imread(name);
            initial_segment(initial_segment < 255) = 0;
        end
    end
    
    % get the initial flows
    [m, n, c] = size(video_frames{1});
    n_frames = size(video_frames,2);
    initial_flows = cell(1, n_frames-1);
    optical_flow = vision.OpticalFlow('OutputValue', ...
        'Horizontal and vertical components in complex form');
    P = step(optical_flow, double(rgb2gray(video_frames{1})));
    for j =2:n_frames
        initial_flows{j-1} = step(optical_flow, ...
            double(rgb2gray(video_frames{j})));
    end
    
    % convert flows to u,v and then to max likely weights
    initial_weights = cell(1,size(initial_flows,2));
    for j = 1:(n_frames-1)
        fprintf('Computing initial weights for frame %d...\n', j);
        
        U = real(initial_flows{j});
        V = imag(initial_flows{j});
        
        % clamp / quantize U and V
        U(U >  b) =  b;
        U(U < -b) = -b;
        V(V >  b) =  b;
        V(V < -b) = -b;
        U = int8(U);
        V = int8(V);
        
        weight_map = zeros(m, n, 2*b+1, 2+b+1);
        for x = 1:m
            for y = 1:n
                weight_map(x,y,U(i,j)+b+1,V(i,j)+b+1) = 1;
            end
        end
        initial_weights{j} = weight_map;
    end
end
