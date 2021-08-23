%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction Code
% Based on Excel from myVGL export file
% The problem still can not indicater the pore size and location
% Avoid thememopry problem
% visulization has been done
% Densily and porosiyt calculation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear
close all
%% read the file
x = xlsread('Pore.xlsx');
%% Parameter setup
Voxel_size = 18.604;
meshsize = 0.018604; % 20um 0.02 mm
% object_dimension (micrometer)
Height = 2000;
Width  = 2000;
Long = 60000;
VolumeObject = Height * Width * Long;
%% creat the data matrix
diameter = x(:,3);
voxel = x(:,8);
volume = x(:,7);
%% realcoordinate
% x_dir = x(:,4);
% y_dir = x(:,5);
% z_dir = x(:,6);
%% the coordinate  z->x, y->z, x->y
x_dir_soft = x(:,14);
y_dir_soft = x(:,15);
z_dir_soft = x(:,16);
% soft_ware z->x, soft_ware y->z, soft_ware x->y
x_dir = z_dir_soft;
y_dir = x_dir_soft;
z_dir = y_dir_soft;

%% projection coordinate
% x_dir = x(:,17);
% y_dir = x(:,18);
% z_dir = x(:,19);

% Data_matrix = [round(x_dir/meshsize),round(y_dir/meshsize),round(z_dir/meshsize),voxel,diameter];z->x, y->z, x->y
Data_matrix = [x_dir,y_dir,z_dir,voxel, diameter,volume];
%% Data_pre-processing
min_dir = [];
for index = 1:3
    min_dir(index) = min(Data_matrix(:,index));
end
for index = 1:3
    for i = 1:size(Data_matrix,1)
        Data_matrix(i,index) = Data_matrix(i,index) - min_dir(index);
    end
end
Data_matrix(Data_matrix(:,5)<Voxel_size,:) = [];

%% create the empty mesh cubic based on the coridnate
% The object dimention 60 mm x 20 mm x 20 mm
% The resolution of x-ray is 18 um 
% the Empty mesh cubic is 3400 point x 1200 point x 1200 point
% The memory problem !!! cnat not use this to show the layer
x = 1:0.01:60;
y = 1:0.01:60; % lim is 20
z = 1:0.01:60; % lim is 20
%% Visualization
figure,scatter3(Data_matrix(:,1),Data_matrix(:,2),Data_matrix(:,3),40,Data_matrix(:,5),'filled');
colorbar('eastoutside')
xlabel('x-direction','Fontsize',40)
ylabel('y-direction','Fontsize',40)
zlabel('z-direction','Fontsize',40)
xlim([2000,23500])
ylim([0,2000])
zlim([0,2000])
set(gca,'Fontsize',20)

% Relative density calculation
Relative_matrix = ((VolumeObject - sum(Data_matrix(:,6)))/VolumeObject)*100;


