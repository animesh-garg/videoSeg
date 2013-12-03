function [] = visualizeFlow(img,U,V)
%VISUALIZEFLOW Given the image and flow fields U,V for each pixel, this
%visualizes the flow field
%   Detailed explanation goes here

interval = 1;

[dimy,dimx,~] = size(img);
x = 1:interval:dimx;
y = 1:interval:dimy;
[X,Y] = meshgrid(x,y);

imagesc(img);
colormap('gray')
hold on;
quiver(X,Y,U(1:interval:dimy,1:interval:dimx),V(1:interval:dimy,1:interval:dimx));

end

