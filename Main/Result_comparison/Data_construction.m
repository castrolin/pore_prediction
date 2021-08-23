%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based the projection image
% Rotating the data point
% create: 2021/6/8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear
close all

image = imread('Scan2\test300114.tif');
image_double = double(image);
% figure,
%     subplot(1,2,1)
%     plot(image_double(:,500))
%     subplot(1,2,2)
%     plot(image_double(100,:))
% for i = 100:100:1000
%     figure(round(i/100)),
%     subplot(1,2,1)
%     plot(image_double(:,i))
%     subplot(1,2,2)
%     plot(image_double(i,:))
% end
s = zeros(1000,1000);
for i = 1:1000
    for j =1:1000
        s(i,j) = round(255* sin(255*(3.14/180)*j)) ;
    end
end
% figure,imshow(uint8(s))
 %% testing rotate image
 s_new = zeros(1000,1000);
 theta = 360;
 for y_new = 1:1000
     for x_new = 1:1000
         i = round(abs(x_new * cos(theta) + y_new * sin(theta))+0.5);
         j = round(abs(-x_new * sin(theta) + y_new * cos(theta))+0.5);
         if i>=0 && i< 1000 && j>=0 && j<1000
             s_new(y_new,x_new) = s(j,i);
         else
             s_new(y_new,x_new) =0;
         end
     end
 end
 
 % figure,imshow(uint8(s_new))
 
 %% simulation the single array( a lot)

 sim = zeros(1000,1000);
 for i = 1:1000
     sim(500,i) = sin(i);
 end
 n = 1;
 sim_new = zeros(1000,1000,25);
 
for theta = 0:15:360
    for y_new = 1:1000
        for x_new = 1:1000
            i = round(abs(x_new * cos(theta) + y_new * sin(theta))+0.5);
            j = round(abs(-x_new * sin(theta) + y_new * cos(theta))+0.5);
            if i >=0 && i<1000 && j>=0 && j<1000
                sim_new(y_new,x_new,n) = sim(j,i);
            else
                sim_new(y_new,x_new,n) = 0;
            end
        end
    end
    n = n+1;
end
                
% for index = 1: size(sim_new,3)
%     figure(index),imshow(sim_new(:,:,index))
% end
sim_add = zeros(1000,1000);
for  index = 1: size(sim_new,3)
    sim_add = sim_add + sim_new(:,:,index);
end
figure,imshow(sim_add);

load tetmesh

figure,scatter3(X(:,1),X(:,2),X(:,3));
figure,trisurf(tet,X(:,1),X(:,2),X(:,3));
%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear,close all

cla
load accidents hwydata                             % load data

long = -hwydata(:,2);                              % longitude data
lat = hwydata(:,3);                                % latitude data
rural = 100 - hwydata(:,17);                       % percent rural data
fatalities = hwydata(:,11);                        % fatalities data

scatter3(long,lat,rural,40,fatalities,'filled')    % draw the scatter plot
ax = gca;
ax.XDir = 'reverse';
view(-31,14)
xlabel('W. Longitude')
ylabel('N. Latitude')
zlabel('% Rural Population')

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Fatalities per 100M vehicle-miles';

%%%%%%%%%%%%%%%%%%%%
% clc,clear,close all
% X = 1:(1/20):60;
% Y = 1:(1/20):60;
% Z = 1:(1/20):60;
% dd = zeros(1,size(X,2));
% for ii = 1:size(X,2)
%     if mod(ii,9) == 0
%         dd(ii) = ii;
%     end
% end
% figure,scatter3(X(:),Y(:),Z(:),40,dd,'filled')
% xlim([0,60])
% ylim([0,20]);
% zlim([0,20]);



