clc;clear;close;
% #2. 2d Image signal filtering

% 2-1. Sketch Magnitude of original image at frequency domain 
img = im2double(imread('cameraman.jpg'));  % read image
IMG = fftshift(fft2(img)); % Fourier Transform of image
figure(1),
subplot(121), imshow(log(1+abs(IMG)), []); title('Magnitude Spectrum of original image'); % frequency domain 
subplot(122), imshow(angle(IMG), []); title('Phase Spectrum of original image');


% 2-2. ILPF in 2D domain
[M,N] = size(img); % The size of the img in pixels
% M : number of rows (height of the img)
% N : number of columns (width of the img)

% ILPF in 2D domain
D0 = 45; % pi/4 Cutoff Frequency (window size:64) <=> radius of the disk
u = 0:(M-1); % u vector 범위 설정
v = 0:(N-1); % N개 v vector 범위 설정
idx = find(u>M/2); % When u is greater than half the number of rows
u(idx) = u(idx)-M; % subtract M
idy = find(v>N/2); % When v is greater than half the number ofcolumns
v(idy) = v(idy)-N; % subtract N
[U,V] = meshgrid(u, v); % 2D grid which contains the coordinates of vectors u and v. 
D = sqrt(U.^2+V.^2); % Euclidean Distance
H_L = fftshift(double(D <= D0)); % Determining the filter mask by comparing with cutoff freq
figure(2),imshow(H_L), title('Ideal LowPass Filter');


% 2-3. IHPF in 2D domain
H_H = 1-H_L; % Transformation in frequency domain 
Y_L = H_L.*IMG; % Ideal Lowpass Filtering
y_L = real(ifft2(ifftshift(Y_L))); % Inverse FT of Y_L
Y_H = H_H.*IMG; % Ideal Highpass Filtering
y_H = real(ifft2(ifftshift(Y_H))); % Inverse FT of Y_H
figure(3),
subplot(121), imshow(y_L, [ ]); title('Output image of ILPF');
subplot(122), imshow(y_H, [ ]); title('Output image of IHPF');
%imwrite(y_L,'out_ILPF.png'); imwrite(y_H,'out_IHPF.png'); 

% 2-4. Gaussian filter and Laplacian filter
% (1) Gaussian filter (LPF)
% Gaussian filter kernel
sigma = 1;
kernal_size = 3;
rad = floor(kernal_size/2);
[a, b] = meshgrid(-rad:rad, -rad:rad);
kernel_G = (1/(2*pi*(sigma.^2)))*exp(-1 * (a.^2 + b.^2)./(2*sigma*sigma));
kernel_G = kernel_G / sum(kernel_G(:));  % Normalization
% zero padding for gaussian kernel
zp_kernel = zeros(size(img)) ;
zp_kernel(1:kernal_size,1:kernal_size) = kernel_G;
ZP_kernel = fftshift(fft2(zp_kernel)); % FT

% (2) Laplacian filter (HPF)
KL = ones(size(img));
kernel_L = [0 -1 0; -1 4 -1; 0 -1 0]; 
% zero padding for laplacian kernel
zp_kernel_L = zeros(size(img)) ;
zp_kernel_L(1:3,1:3) = kernel_L;
ZP_kernel_L = fftshift(fft2(zp_kernel_L));

figure(4),
subplot(121),imshow(log(1+abs(ZP_kernel)), [ ]); title('Gaussian filter');
subplot(122),imshow(log(1+abs(ZP_kernel_L)),[ ]); title('Laplacian filter');


% 2-5. Filtering the noised image by ILPF / GLPF
noised_img = im2double(imread('noised_img.jpg.png')); 
NOISED_IMG = fftshift(fft2(noised_img));  % Fourier Transform of noised_img

filtered_ILPF = H_L.*NOISED_IMG; % ILPF Filtering in freq domain
filtered_ilfp = real(ifft2(ifftshift(filtered_ILPF))); % Inverse FT of filtered_ILPF
filtered_gaussian = conv2(noised_img, kernel_G);  % Gaussian Filtering in image domain

figure(5),
subplot(121), imshow(filtered_gaussian, [ ]); title('Output image of Gaussian LPF (noised)');
subplot(122), imshow(filtered_ilfp, [ ]); title('Output image of ILPF (noised)');
%imwrite(filtered_gaussian,'out_Gaussian.png'); imwrite(filtered_ilfp,'out_ILPF_noised.png');

% Comparison the performance of LPF using PSNR metric
img_g = imresize(img,[258,258]);
psnr_ILPF = psnr(filtered_gaussian, img_g);
psnr_gaussian = psnr(filtered_ilfp, img);
display(psnr_gaussian); 
display(psnr_ILPF);


% 2-6. Filtering the original image on frequency domain(multiplication) by IHPF / LHPF
Y_Lap = ZP_kernel_L.*IMG; 
y_Lap = real(ifft2(ifftshift(Y_Lap)));

figure(6),
subplot(121), imshow(y_Lap, [ ]); title('Output image of Laplacian HPF');
subplot(122), imshow(y_H, [ ]); title('Output image of IHPF');
%imwrite(y_Lap,'out_Laplacian.png');


