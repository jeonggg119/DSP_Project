% 1-1. x[n] and |X(w)| for 'x1.wav'
clc;clear;  

[x,F] = audioread('x1.wav'); 
n = length(x);  % length of signal 
w = linspace(-pi,pi,n);  % set frequency domain interval
X = fftshift(fft(x));  % FT of signal and make freq 0 be highest magnitude (Move center)

figure(1),
subplot(121), plot(x); xlabel('n'); ylabel('x[n]'); title('Original signal x[n]'); ylim([-1.2 1.2]);
subplot(122), plot(w,log(1+abs(X))); xlabel('w'); ylabel('|X(w)|'); title('Magnitude of original signal(|X(w)|)');

% 1-2. Ideal Low-pass filter 
H_L = 1.*(w<=pi/80 & w>=-pi/80) + 0.*(w>pi/80 & w<-pi/80);  % ILPF (저주파 영역만 1로 설정, 나머지 영역은 0으로 설정) Frequency response
figure(2),
subplot(211), plot(w,abs(H_L)); xlabel('w'); ylabel('|H(w)|'); title('Magnitude of H(w)(IDLF)'); ylim([0 0.75]);
subplot(212), plot(w,angle(H_L)); xlabel('w'); ylabel('phase of H(w)'); title('Phase of H(w)(IDLF)');

% 1-3. IDLF Filtering
Y_L = X.*H_L';  % Filtering = LTI system
y_L = real(ifft(ifftshift(Y_L)));  % Inverse FT of output signal and move center of it, real for 부동 소수점 유효숫자 오차 제거
figure(3),
subplot(121), plot(y_L); xlabel('n'); ylabel('y[n]'); title('Output signal y[n](IDLF)'); ylim([-1.2 1.2]);
subplot(122), plot(w,log(1+abs(Y_L))); xlabel('w'); ylabel('|Y(w)|'); title('Magnitude of output signal(IDLF)'); 
% audiowrite('lowest.wav',y_L,F);
% sound(y_L,F);

% 1-4. Ideal Band-pass filter 
H_B = 1.*( w>=-pi/40 & w<=-pi/65) + 0.*(w>-pi/65 & w<pi/65) + 1.*(w>=pi/65 & w<=pi/40); % IBPF (원하는 middle 영역만 0을 기준으로 대칭적 추출) Frequency response
figure(4),
subplot(211), plot(w,abs(H_B)); xlabel('w'); ylabel('|H(w)|'); title('Magnitude of H(w)(IDBF)'); ylim([0 0.75]);
subplot(212), plot(w,angle(H_B)); xlabel('w'); ylabel('phase of H(w)'); title('Phase of H(w)(IDBF)');

% 1-5. IDBF Filtering
Y_B = X.*H_B';  % Filtering = LTI system 
y_B = real(ifft(ifftshift(Y_B)));  % Inverse FT of output signal, real for 부동 소수점 유효숫자 오차 제거
figure(5),
subplot(121), plot(y_B); xlabel('n'); ylabel('y[n]'); title('Output signal y[n](IDBF)'); ylim([-1.2 1.2]);

subplot(122), plot(w,log(1+abs(Y_B))); xlabel('w'); ylabel('|Y(w)|'); title('Magnitude of output signal(IDBF)');
% audiowrite('middle.wav',y_B,F); 
% sound(y_B,F);