% %% plot length
close all
figure,
length=lenavg;
error=lenstd;
h=barwitherr(error,length);
set(gca,'XTickLabel',{'1','2','3','4','5','6','7'},'FontSize',25)
%legend('Monitoring',30)
xlabel('Clustering Number','FontSize',25)
ylabel('\mum','FontSize',25)
ylim([0 500]) % specific y axis range
title('extracted length for scaning speed 1000mm/s','FontSize',25)

%% plot width
figure;
width=widavg;
error=widstd;%error bar
h=barwitherr(error,width);
set(gca,'XTickLabel',{'1','2','3','4','5','6','7'},'FontSize',25)
%legend('experiment','Simon simulation','Eason')
xlabel('Clustering Number','FontSize',25)
ylabel('\mum','FontSize',25)
ylim([0 500]) % specific y axis range
title('extracted width for scaning speed 1000mm/s','FontSize',25)

%% plot ratio
figure;
ratio=ratioavg;
% error=ratiostd;%error bar
error = [0,0,0,0,0,0,0];
h=barwitherr(error,ratio);
set(gca,'XTickLabel',{'1','2','3','4','5','6','7'},'FontSize',25)
%legend('experiment','Simon simulation','Eason')
xlabel('Clustering Number','FontSize',25)
ylabel('dimensionless')
title('extracted ratio for scaning speed 1000mm/s')
%% plot spatter
% figure;
% width=[49 42 24 23;15 16 12 6;12 8 4 5;7 6 4 4; 5 4 8 2];
% error=zeros(size(width));
% h=barwitherr(error,width)
% set(gca,'XTickLabel',{'0~50','50~75um','75~100um','100~125um','>125um'})
% legend('300W','250W','200W','150W');
% xlabel('spot diameter')
% ylabel('count')
% title('histogram for different size')