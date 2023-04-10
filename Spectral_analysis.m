% Spectral Analysis of an image 
%Turn off terrain for Ortho-image (flat image) 

image = imread( 'White_Sands_GE_10_2009.png');
whos image %Determine size of image - Will provide size & bytes 
imshow(image)

%% Display image in black and white

image_rgb = rgb2gray(image); 
imshow(image_rgb)

%% Define image parameters 
[x,y] = ginput(2);  %Make a live image

%204 pixels = 481 m

% 204/481 m  = 0.42 m ; F_n = 0.42/2 = 0.21 m (nyquist) -> lambda = 1/.21 m
% = 5 m - smallest spatial resolution

%Use imagesc for nice color map and scale on axis

%% Create axis with 2.4 m per pixel

image = imread('White_Sands_GE_10_2009.png');
resolution = 2.4; 
xrange = size(image,2) * resolution; % in meters
yrange = size(image,1) * resolution; % in meters
axis_x = 0:resolution:xrange;
axis_y = 0:resolution:yrange;

imagesc(axis_x, axis_y, image);
xlabel('Distance (m)');
ylabel('Distance (m)');

%%
subplot(2,1,1)
data = image(200,:); 
plot(data,'m'); 

subplot(2,1,2)
data1 = image(300, :); 
plot(data1, 'b')

%% Create an annotation tool
 image = imread( 'White_Sands_GE_10_2009.png');
whos image %Determine size of image - Will provide size & bytes 
imagesc(image)

[x,y] = ginput(2); 
plot(x,y, 'mo'); 
n = 51; 
xps = linspace(x(1), x(2), n); 
yps = linspace(y(1), y(2), n)
plot(xps, yps, '-b')

%% Interpolate

zps = interp2(xs, ys, image, xps, yps); 
plot(zps)