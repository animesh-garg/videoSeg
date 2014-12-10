function [init_image] = read_image(path)
init_image = imread (path);
init_image = logical  (init_image);
end