clc;clear;close;
% 2-1. Sketch Magnitude and Phase of original image at frequency domain 
img = im2double(imread('twigs.png'));  % read image
X = fftshift(fft2(img));  % Fourier Transform of image

figure(1),
subplot(121), imshow(log(1+abs(X)), []); title('Magnitude Spectrum of original image'); % frequency domain 
% middle white area : highest magnitude = image의 DC 성분 = F(0,0), at low freq.
% 동서남북 방향 bright area : noise -> 제거할 부분
subplot(122), imshow(angle(X), []); title('Phase Spectrum of original image');


% 2-2. Design Filter to remove noise in the image and Reconstruct the image
dcm = datacursormode; % using datacursor to find the location of noises.
dcm.Enable = 'on';
dcm.DisplayStyle = 'window';

% filter of clearing noise by multiplying noise area to zero point
H = ones(256,256); % set filter size same as image size
H(1:1:88, 128:1:131)=0; % clear up_side noise 
H(168:1:256, 128:1:132)=0; % clear down_side noise 
H(126:1:131, 1:1:88)=0; % clear left_side noise
H(126:1:132, 170:1:256)=0; % clear right_side noise

figure(2),
subplot(121), imshow(abs(H)); title('Magnitude Spectrum of Filter'); % frequency domain
subplot(122), imshow(angle(H)); title('Phase Spectrum of Filter');
  
Y = X.*H;  % LTI system Filtering
y = real(ifft2(fftshift(Y)));  % Inverse FT of output signal(denoised image)


% 2-3. Sketch Magnitude and Phase of denoised image at frequency domain 
figure(3),
subplot(121), imshow(log(1+abs(Y)),[]); title('Magnitude Spectrum of denoised image');
subplot(122), imshow(angle(Y),[]); title('Phase Spectrum of denoised image');

figure(4), % compare images to check filtering
subplot(121), imshow(img); title('Original image');
subplot(122), imshow(y); title('Denoised image');
imwrite(y,'denoisedimage.png');
