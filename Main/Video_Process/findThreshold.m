function [T,boundary]=findThreshold(I,boundary)
pixel=17;
horizontalProjection=sum(I');
centerY=find(horizontalProjection==max(horizontalProjection));
centerY=centerY(1);
y_raw=I(centerY,:);
y_raw=double(y_raw);
y=fourierLowPass(y_raw,30,10);
d1=diff(y);
d2=diff(d1);
d2min=find(d2==min(d2));
x=1*pixel:pixel:256*pixel;
spline=fit(x',y_raw','smoothingspline','SmoothingParam',1.3e-6);
y_spline=spline([1*pixel:pixel:256*pixel]);
TF = islocalmin(d2(1:d2min));

%% plot raw, 1st, 2nd, low pass filter
% figure();
% plot(1*pixel:pixel:256*pixel,double(y_raw),'linewidth',5);
% hold on
% plot(1*pixel:pixel:256*pixel,y,'linewidth',5);
d1=diff(y);
d2=diff(d1);
% plot(1*pixel:pixel:255*pixel,d1*10,'linewidth',5);
% plot(1*pixel:pixel:254*pixel,d2*100,'linewidth',5);
% grid on
% xlabel('length[um]','FontSize', 30)
% ylabel('intensity','FontSize', 30)
% legend({'raw data','after low pass filter','first derivative','second derivative'},'FontSize',14);
% set(gca,'FontSize',30)
% set(gca,'linewidth',3)
%% plot raw, 2nd, spline fitting
% figure();
% plot(1*pixel:pixel:256*pixel,double(y_raw),'linewidth',5);
% hold on
% plot(1*pixel:pixel:256*pixel,y,'linewidth',5);
% plot(1*pixel:pixel:256*pixel,y_spline,'linewidth',5);
% legend1=legend({'raw data','low pass filter','spline fitting'},'FontSize',30);
% set(legend1,'Location','northwest','FontSize',30);
% grid on;
% hold off;
% xlabel('length[um]','FontSize', 30)
% ylabel('intensity','FontSize', 30)
% set(gca,'FontSize',30)
% figure()
% d1_spline=diff(y_spline);
% d2_spline=diff(d1_spline);
% plot(1*pixel:pixel:254*pixel,d2*100,'linewidth',5);
% hold on
% plot(1*pixel:pixel:254*pixel,d2_spline*100,'linewidth',5);
% grid on
% hold off;
% xlabel('length[um]','FontSize', 30)
% ylabel('intensity','FontSize', 30)
% legend2=legend({'2nd derivative for low pass filter','2nd derivative for smooting spline'},'FontSize',30);
% set(legend2,'Location','northwest','FontSize',30);
% set(gca,'FontSize',30)
%% 將monitoring結果與simulation 重疊

p='E:\codenew\250-600.csv';
sim=xlsread(p);
sim=sim(:,4);
sim=smooth(sim);
peak=find(y==max(y));
% try
% 150W  a=90  b=+17 p=1.6 %p為output到excel中一個資料點代表的距離 ex 1000模擬點表示1.6mm-->1point=1.6um 
% 200W  a=80  b=+27 p=1.6
% 250W  a=81  b=26  p=1.9
% 300W  a=74  b=33 p=1.9
try
a=70;  %a+b = 107
b=37;
p=1.7; 
yy=y(peak-a:peak+b);
yyy=resample(yy,1000,108);  %resample(input matrix, new scale length, original scale length)
d11=d1(peak-a:peak+b);
d11=resample(d11,999,108);
d21=d2(peak-a:peak+b);
d21=resample(d21,998,108);
d1u=diff(sim);
d2u=diff(d1u);

figure(1);
yyaxis left
plot(1*p:p:1000*p,yyy,'linewidth',5)
ylabel('intensity','FontSize', 30)
yyaxis right
plot(1*p:p:1000*p,sim,'linewidth',5);
ylabel('temperature[K]','FontSize', 30)
xlabel('position[um]','FontSize', 30)
legend1=legend({'monitored intensity','simulated temperature'},'FontSize',30);
set(legend1,'Location','northwest','FontSize',30);
title('overlap of the temperature and intensity profile ')
grid on;
set(gca,'FontSize',30)

% figure(2);
% yyaxis left
% plot(1*1.6:1.6:999*1.6,d11)
% yyaxis right
% plot(1*1.6:1.6:999*1.6,d1u);
% legend('experiment','simulation')
% xlabel('position');
% title('一階微分')

figure(3)
yyaxis left
plot(1*p:p:998*p,d21,'linewidth',5)
ylabel('arbitrary unit','FontSize', 30)
yyaxis right
plot(1*p:p:998*p,smooth(d2u),'linewidth',5);
ylabel('arbitrary unit','FontSize', 30)
legend1=legend('monitored data','simulated data');
set(legend1,'Location','northwest','FontSize',30);
xlabel('position[um]','FontSize', 30)
title('overlap of the second derivative data ')
grid on;
set(gca,'FontSize',30)

end
% 
% 
% 
localmin=find(TF==1);
if nargin==1
    boundary=localmin(end);
    T=y(boundary)/255;
    close all;
else

    [~, index] = min(abs(localmin-boundary));
    if boundary>localmin(index)
        T=y(boundary)/255;
        close all;
    else
        boundary=localmin(index);
        T=y(boundary)/255;
        close all;
    end

end


