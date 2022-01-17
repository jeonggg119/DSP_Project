clc;clear;close;
% #1. Audio signal filtering

% 1-1. x[n] and |X(w)| for 'x[n].wav' &  Generate White Gaussian noise v[n]
[x,F] = audioread('x[n].wav'); 
w = linspace(-pi,pi,length(x));  % set frequency domain interval
X = fftshift(fft(x));  % FT of signal and make freq 0 be highest magnitude (Move center)

figure(1),
subplot(121), plot(x); xlabel('n'); ylabel('x[n]'); title('Original voice signal x[n]');
subplot(122), plot(w,log(1+abs(X))); xlabel('w'); ylabel('|X(w)|'); title('Magnitude of original voice signal');
v = 0.02 * randn(size(x)); % White Gaussian noise N(0,0.02)
% sound(v+x,F)


% 1-2. Truncated ideal impulse response h_d at time domain
N =39;
r = length(x)/2;
n = -r:1:r-1;
h_d = ((0.25*sinc(0.25*n)).*(n>=0 & n<=38) + 0.*(n<0 & n>38))';
H_d = fftshift(fft(h_d));

figure(2),
plot(w,real(H_d)); ylim([0 0.7])
xlabel('w'); ylabel('real(H_d)'); title('Gibbs phenomenon by Truncated IIR');
% Gibbs phenomenon : oscillation near discountinuites


% 1-3. FIR lowpass filter h[n] with 39 point Hanning Window
h = ((0.5*sinc(0.5*n)).*(n>=0 & n<=38) + 0.*(n<0 & n>38))';  % using normalized frequency
n_shift = (39-1)/2; % (N-1)/2 linear phase time shift
window = (0.5.*(1-cos(2*pi*(n-n_shift)/38)).*(n>=0 & n<=38) + 0.*(n<0 & n>38))';  % shifted window

h_fir = h.*window ; % Time domain multiplication -> FIR filter
H_fir = fftshift(fft(h_fir));  % FT of h_fir

figure(3),
plot(w,log(1+abs(H_fir))); xlabel('w'); ylabel('|H(w)|'); title('Magnitude of FIR lowpass filter');
% peak sidse lobe amplitude 감소(good), approx width of main lobe 증가(bad)


% 1-4. Filtering v[n] with H(w)
V = fftshift(fft(v)); % Fourier transform of v[n]
Vf = V.*H_fir; % Filtering noise 
vf = real(ifft(ifftshift(Vf))); % Inverse FT of Vf(w)

figure(4),
subplot(121), plot(vf); xlabel('n'); ylabel('vf[n]'); title('filtered noise signal vf[n]'); 
subplot(122), plot(w,log(1+abs(Vf))); xlabel('w'); ylabel('|Vf(w)|'); title('Magnitude of filtered noise signal'); 

x_d = x + vf; % Original sound + filtered noise sound by hanning window method 
X_d = fftshift(fft(x_d)); 
% sound(x_d,F);


% 1-5. My own filter H2(w) using pole-zero placement ->  전체 신호를 필터링 (X_d)
% Where I want to Emphasize (using poles) : z = 0.5, -0.5
% Where I want to Deemphasize (using zeros) : 0 (double)
% Poles : always inside the unit circle

z = exp(j*w);
H_W_2=( 0.765*  ((z).^2)./( (z-0.45).*(z+0.45) ) )'; % My own filter

X_d2 = X_d.*H_W_2;  % Filtering
x_d2 = real(ifft(ifftshift(X_d2)));  % Inverse FT of X_d2
% sound(x_d2,F);
% audiowrite('out_myfilter.wav',x_d2,F);

figure(5),
subplot(121),plot(w,log(1+abs(H_W_2))); xlabel('w'); ylabel('|H_2(w)|'); title('Magnitude of My own filter');
subplot(122),plot(w,angle(H_W_2)); xlabel('w'); ylabel('phase of H_2(w)'); title('Phase of My own filter');

% The error score of given equation
sum = 0;
for i = 1 : length(x)-1
  sum = sum + (x(i)-x_d2(i)).^2; 
end
score = sqrt(sum);
display(score);




% figure(6), %(my own) for comparing orginal signal and filtered noise signal
% plot(w,log(1+abs(X))); 
% hold on;
% plot(w,log(1+abs(X_d))); xlabel('Vf'); ylabel('|Vf(w)|'); title('Magnitude of filtered noise signal');
% title('Magnitude of orginal signal and filtered noise signal by window method')
% xlabel('w'); ylabel('Magnitude'); legend('original, |X(w)|','filtered, |Vf(w)|');
% 
% figure(7), %(my own) for checking output signal of my filter
% subplot(121), plot(x_d2);  xlabel('n'); ylabel('x_d2[n]'); title('Output signal of my filter');
% subplot(122), plot(w,log(1+abs(X_d2))); xlabel('w'); ylabel('|X(w)|'); title('Magnitude of Output signal of my filter');
% 
% figure(8), %(my own) for comparing orginal signal and filtered noise signal by my filter
% plot(w,log(1+abs(X))); xlabel('w'); ylabel('|X(w)|'); title('Magnitude of original voice signal');
% hold on;
% plot(w,log(1+abs(X_d2))); xlabel('w'); ylabel('|X_d2(w)|'); title('Magnitude of filtered noise signal');
% title('Magnitude of orginal signal and filtered noise signal by my filter')
% xlabel('w'); ylabel('Magnitude'); legend('original, |X(w)|','filtered, |X_d2(w)|');





