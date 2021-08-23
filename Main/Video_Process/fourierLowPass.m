% 
% t = linspace(1,10,1024); 
% x = -(t-5).^2 + 2; 
% y = awgn(x,0.5); 
% Y = fft(y,1024); 
% 
% r = 20; % range of frequencies we want to preserve 
% 
% rectangle = zeros(size(Y)); 
% rectangle(1:r+1) = 1;    % preserve low +ve frequencies 
% y_half = ifft(Y.*rectangle,1024); % +ve low-pass filtered signal 
% rectangle(end-r+1:end) = 1;   % preserve low -ve frequencies 
% y_rect = ifft(Y.*rectangle,1024); % full low-pass filtered signal 
% 
% hold on; 
% plot(t,y,'g--');
% plot(t,x,'k','LineWidth',2);
% plot(t,y_half,'b','LineWidth',2);
% plot(t,y_rect,'r','LineWidth',2); 
% legend('noisy signal','true signal','+ve low-pass','full low-pass','Location','southwest') 
% % clc,clear;
% % cdata=imread('D:\Works\SLM melt pool\code\testing image\snapshot.bmp');
function y_gauss=fourierLowPass(y,r,sigma)

Y=fft(y);

s = size(y,2);      %¯x°}ªø«×¤@­P

gauss = zeros(size(Y)); 
gauss(1:r+1) = exp(-(1:r+1).^ 2/(2 * sigma^2)); % +ve frequencies 
gauss(end-r+1:end) = fliplr(gauss(2:r+1));   % -ve frequencies 
y_gauss = ifft(Y.*gauss,s); 

