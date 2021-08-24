%% Analysis the cluster data based on likelihood 
% Length Width Ratio spatter distribution 
% Data was extracted from excel file the  title is Test_(name of the layer).xlsx
% The memory issu has to be noticed

% Note:
% data{1,1} = centroid of cluster info
% data{1,2} = Number of each cluster
% data{1,3} = mean and std for [Length, width, Ratio, Angle, Number of cluster]
% data{1,4} = index of melt pool in each cluster

clc,clear
close all
%% read excel and post-processing for cluster data
path = 'X:\Castro\statistical_analysis\ecel_0715\Test100.xlsx';
[~,sheet_name] = xlsfinfo(path);
for k = 1:numel(sheet_name)
    data_sheet{k} = xlsread(path,sheet_name{k});
end
%---------------------------------------------------------------------------
%% data_cleaning
n = 0;
for i = 1:size(data_sheet{1,2},2)
    n = n+1;
end
data_sheet{1,3}(n+1:end,:) = [];
data_sheet{1,1}(find(data_sheet{1,1}(:,1) == 0),:) = [];

%% likelihood function(maximum likelihood function)
Length_mean = [data_sheet{1,3}(:,1)]*2;
Width_mean = [data_sheet{1,3}(:,2)]*2;
Ratio_mean = [data_sheet{1,3}(:,3)];
Angle_mean = [data_sheet{1,1}(:,4)];

% the number of melt pool in each cluster
NOD = data_sheet{1,2};
ns = sum(NOD);
numbins = 10;
%-------------------------------------------------------------------
% Length
figure
for number = 1: size(NOD,2)
    mu = Length_mean(number); %Population parameter
    n = NOD(number);
    rng('default');
    samples = exprnd(mu,n,ns);
    means = mean(samples);
    [phat,pci] = mle(means); % Phat(1) = mean, Phat(2) = std
    histogram(means,numbins,'Normalization','pdf');
    hold on
    x = min(means):max(means);
    y = normpdf(x,phat(1),phat(2));
    plot(x,y,'g','LineWidth',2)
end
xlabel('Length(\mum)')
ylabel('PDF')
set(gca,'FontSize',20)
%-------------------------------------------------------------
% Width
figure
for number = 1: size(NOD,2)
    mu = Width_mean(number); %Population parameter
    n = NOD(number);
    rng('default');
    samples = exprnd(mu,n,ns);
    means = mean(samples);
    [phat,pci] = mle(means); % Phat(1) = mean, Phat(2) = std
    histogram(means,numbins,'Normalization','pdf');
    hold on
    x = min(means):max(means);
    y = normpdf(x,phat(1),phat(2));
    plot(x,y,'g','LineWidth',2)
end
xlabel('Width(\mum)')
ylabel('PDF')
set(gca,'FontSize',20)

%---------------------------------------------------------------
% Ratio
figure
for number = 1: size(NOD,2)
    mu = Ratio_mean(number); %Population parameter
    n = NOD(number);
    rng('default');
    samples = exprnd(mu,n,ns);
    means = mean(samples);
    [phat,pci] = mle(means); % Phat(1) = mean, Phat(2) = std
    histogram(means,numbins,'Normalization','pdf');
    hold on
    x = min(means):0.001:max(means);
    y = normpdf(x,phat(1),phat(2));
    plot(x,y,'g','LineWidth',2)
end
xlabel('Ratio(W/L)')
ylabel('PDF')
set(gca,'FontSize',20)

%------------------------------------------------------------------
% Angle 
% the number of spatter in each cluster
NOS = data_sheet{1,3}(:,5);

Each_Number = round(NOS + data_sheet{1,3}(:,10));

figure
for number = 1: size(NOS,1)
    mu = Angle_mean(number); %Population parameter
    n = Each_Number(number);
    rng('default');
    samples = exprnd(mu,n,sum(Each_Number));
    means = mean(samples);
    %[phat,pci] = mle(means,'distribution','tLocationScale'); % Phat(1) = mean, Phat(2) = std
    histogram(mu,size(NOS,1),'Normalization','pdf');
    %grid on;
    xlim([0, 30]);
    hold on
%     x = 0:0.05:(mu + data_sheet{1,3}(number,10));
%     y = normpdf(x,mu,data_sheet{1,3}(number,10));
%     plot(x,y,'g','LineWidth',2)
end
xlabel('Angle(\theta)')
ylabel('PDF')
ylim([0,10])
set(gca,'FontSize',20)

disp(['num_spatter:' num2str(round(sum(Each_Number)))])
