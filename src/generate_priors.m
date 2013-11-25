% Script to generate priors for our video segementatino algorithm
%   Uses the moseg_dataset and 3rd party NCuts code
%   Author: Jeff 
function [initial_segment, initial_flows, initial_weights, avg_I_f, avg_I_b] = ...
    generate_priors(window_size, video_folder, downsample_rate)

%% Process all images in the given video sequence
d = downsample_rate; % downsampling values
b = window_size;

% folder files 
folder_files = dir(video_folder);
video_frames = {};
   
% read in all of the frames
k = 1;
for j = 1:size(folder_files,1)
    name = folder_files(j).name;
    if length(name) > 4 && ...
        strcmp(name((size(name,2)-3):size(name,2)), '.jpg')
        
        name = sprintf('%s/%s',video_folder, name);
        fprintf('Loading image %s ...\n',name);
        video_frames{k} = imread(name);
        k = k+1;
    elseif length(name) > 6 && ...
        strcmp(name((size(name,2)-5):size(name,2)), '01.pgm')
        name = sprintf('%s/%s',video_folder, name);
        initial_segment = imread(name);
        
        % pick out only one of the possible object to segment, and store
        % initial segmentation as 0 or 1 valued array
        initial_segment(initial_segment < 255) = 0;      
        initial_segment(initial_segment==255) = 1;
    end
end
    
% get the initial flows
[m, n, c] = size(video_frames{1});
n_frames = size(video_frames,2);
initial_flows = cell(1, n_frames-1);
optical_flow = vision.OpticalFlow('OutputValue', ...
        'Horizontal and vertical components in complex form', ...
        'ReferenceFrameSource', 'Input port', ...
        'Method', 'Lucas-Kanade');
for j =2:n_frames
    initial_flows{j-1} = step(optical_flow, ...
        double(rgb2gray(video_frames{j})), ...
        double(rgb2gray(video_frames{j-1})));
end
    
% convert flows to u,v and then to max likely weights
initial_weights = cell(1,size(initial_flows,2));
for j = 1:(n_frames-1)
    fprintf('Computing initial weights for frame %d...\n', j);
        
    U = real(initial_flows{j});
    V = imag(initial_flows{j});    
    initial_weights{j} = uv_to_weights(U, V, b);
end

%% Compute the appearance model constants
avg_I_f = zeros(1,3);
avg_I_b = zeros(1,3);

[m, n] = size(video_frames{1});
n_frames = size(video_frames,2);
n_pixels_frame = m * n;
n_pixels_foreground = sum(sum(initial_segment));
n_pixels_background = n_pixels_frame - n_pixels_foreground;
n_pixels_total= n_frames * n_pixels_frame;

% get averages for appearance model cost
avg_I_f(1) = (1.0 / n_pixels_foreground) * ...
    sum(sum(video_frames{1}(:,:,1).*initial_segment));
avg_I_f(2) = (1.0 / n_pixels_foreground) * ...
    sum(sum(video_frames{1}(:,:,2).*initial_segment));
avg_I_f(3) = (1.0 / n_pixels_foreground) * ...
    sum(sum(video_frames{1}(:,:,3).*initial_segment));

avg_I_b(1) = (1.0 / n_pixels_background) * ...
    sum(sum(video_frames{1}(:,:,1).*(1-initial_segment)));
avg_I_b(2) = (1.0 / n_pixels_background) * ...
    sum(sum(video_frames{1}(:,:,2).*(1-initial_segment)));
avg_I_b(3) = (1.0 / n_pixels_background) * ...
    sum(sum(video_frames{1}(:,:,3).*(1-initial_segment)));

end
