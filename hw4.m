clc
clear all
close all

image = double(rgb2gray(imread('recorder.jpg')));
imageFFT = fft2(image);
[a,b] = size(imageFFT);

disp('----------Q1 10 Percent-------------------')
compress = 0.1;
[ffta, fftb] = sort(imageFFT(:),'descend');
i=round(a*b*compress);
imagefftcompress = zeros(size(imageFFT));
imagefftcompress(fftb(1:i)) = ffta(1:i);
imagecompress = ifft2(imagefftcompress);
errorfft = norm(imageFFT - imagefftcompress);
errorimage = norm(image - imagecompress);
fprintf('For 10 Percent, L2 norm of error is %.3f, L2 norm of FFT is %.3f \n',errorimage,errorfft)


disp('----------Q1 1 Percent-------------------')
compress = 0.01;
[ffta, fftb] = sort(imageFFT(:),'descend');
i=round(a*b*compress);
imagefftcompress = zeros(size(imageFFT));
imagefftcompress(fftb(1:i)) = ffta(1:i);
imagecompress = ifft2(imagefftcompress);
errorfft = norm(imageFFT - imagefftcompress);
errorimage = norm(image - imagecompress);
fprintf('For 1 Percent, L2 norm of error is %.3f, L2 norm of FFT is %.3f \n',errorimage,errorfft)


disp('----------Q2 a----------------------------')
load r2112.mat;
sound(rush, FS);

figure(1)
i = length(rush);
y = rush;
fs = FS;
d = (-i/2:i/2-1)/i*fs;
value = 2/i*fft(y);
mag = abs(fftshift(value));
plot(d,mag)
title('Q2 a')
xlabel('Freq')
ylabel('Magnitude')
axis([-1000 25000 0 0.03])

disp('----------Q2 b----------------------------')
figure(2)
PSD = mag.*conj(mag);
plot(d, PSD)
title('Q2 b FFT')
xlabel('Freq')
ylabel('Power')
axis([-1000 25000 0 0.0007])
figure(3)
spectrogram(rush,256,120,128,fs)
set(gca,'CLim',[-160,-20]);
title('Q2 b Spectrum')


disp('----------Q2 c----------------------------')
load r2112noisy.mat
i = length(rushnoisy);
noisy = 2/i*fft(rushnoisy);
n = (-i/2:i/2-1)/i*FS;

cleaned = zeros(size(noisy));
range = 1:(i/4);
cleaned(range) = noisy(range)*2;
rushclean = real(ifft(cleaned*i/2));
sound(rush, FS);
sound(rushnoisy, FS);
sound(rushclean, FS);

figure(4);
Pnoisy = abs(noisy.*conj(noisy));
plot(n,fftshift(Pnoisy));
xlabel('Freq');
ylabel('Power');
axis([0 30000 0 1.2]);
title('Q2 c noisy')
figure(5);
Pcleaned = fftshift(abs(cleaned.*conj(cleaned)));
plot(n, Pcleaned)
xlabel('Freq');
ylabel('Power');
axis([-1000 30000 0 0.003]);
title('Q2 c cleaned')

figure(6);
spectrogram(rushnoisy,256,120,128,FS);
set(gca,'CLim',[-160,-20]);
title('Q2 c noisy')
figure(7)
spectrogram(rushclean,256,120,128,FS);
set(gca,'CLim',[-160,-20]);
title('Q2 c cleaned')












