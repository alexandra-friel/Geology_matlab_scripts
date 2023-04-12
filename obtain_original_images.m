

%% Use FFT2 and IFFT2 to obtain original images 
load('slive.mat', 'h', 'y');

H = fft2(h);
Y = fft2(y);

% Compute the Fourier transform of the original image
X = Y ./ H;

% Take the inverse Fourier transform to obtain the original image
x = ifft2(X);

% Display the original and distorted images side by side
sgtitle('FFT2 and inverse FFT to obtain original images')
subplot(1,2,1), imshow(abs(x), []), title('Original Image of Elvis Presley'); 
subplot(1,2,2), imshow(abs(y), []), title('Distorted Image of Elvis Presley');

colormap(pink)

