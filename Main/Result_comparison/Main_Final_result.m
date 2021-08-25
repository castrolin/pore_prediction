%% the main 
% combine read CT result and comparision file 
% Step1. read the exel file and plot the pore size and location
% Step2. read the prediction result and plot the pore size and location
% Step3. compare the prediction with CT
% Calculate the accuracy of model and densiy of CT and prediction

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CT part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read the file
x = xlsread('X:\Castro\recon_pore\excel\Pore_fist_July15.xlsx'); % import CT result(excel file) here
%% Parameter setup
Voxel_size = 18.604;
meshsize = 0.018604; % 20um 0.02 mm
% object_dimension (micrometer) As built dimension
Height = 2000;
Width  = 3000;
Long = 4000;
VolumeObject = Height * Width * Long;

%% create the data matrix
diameter = x(:,3);
voxel = x(:,8);
volume = x(:,7);

%% the coordinate (the excel file -> the title_name: centerx y,z)
% note:The coordinate sytem is different between CT and prediction
y_dir = x(:,6); % soft -> software
z_dir = x(:,4);
x_dir = x(:,5);

Data_matrix = [x_dir,y_dir,z_dir,voxel, diameter,volume];

%% Data_pre-processing (shift the coordinate)
min_dir = [];
for index = 1:3
    min_dir(index) = min(Data_matrix(:,index));
end
for index = 1:3
    for i = 1:size(Data_matrix,1)
        Data_matrix(i,index) = Data_matrix(i,index) - min_dir(index);
    end
end
Data_matrix(Data_matrix(:,4)< 5,:) = [];
figure,scatter3(Data_matrix(:,1),Data_matrix(:,2),Data_matrix(:,3),40,Data_matrix(:,5),'filled');
colorbar('eastoutside')
title('Original datafromCT','Fontsize',40)
xlabel('x-direction','Fontsize',40)
ylabel('y-direction','Fontsize',40)
zlabel('z-direction','Fontsize',40)
% xlim([501,4000])
% ylim([0,2000])
% zlim([0,2000])
set(gca,'Fontsize',20)
% To conduct the FOV issue the start point is from 500 um (x-axis and y-axis)
Data_matrix(Data_matrix(:,1)<1501,:) = [];  % the x-axis position
Data_matrix(Data_matrix(:,2)<2950,:) = []; % the y axis position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the re-scale agian
min_dir = [];
for index = 1:3
    min_dir(index) = min(Data_matrix(:,index));
end
for index = 1:3
    for i = 1:size(Data_matrix,1)
        Data_matrix(i,index) = Data_matrix(i,index) - min_dir(index);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dimenstionto mesh 
Data_matrix_mesh = sortrows([round(Data_matrix(:,1)/Voxel_size),round(Data_matrix(:,2)/Voxel_size),round(Data_matrix(:,3)/35)+1,Data_matrix(:,5)],3);
%% create the empty mesh cubic based on the coridnate
%% Visualization
figure,scatter3(Data_matrix(:,1),Data_matrix(:,2),Data_matrix(:,3),40,Data_matrix(:,5),'filled');
colorbar('eastoutside')
title('CT_inFOV','Fontsize',40)
xlabel('x-direction','Fontsize',40)
ylabel('y-direction','Fontsize',40)
zlabel('z-direction','Fontsize',40)
% xlim([501,4000])
% ylim([0,2000])
% zlim([0,2000])
set(gca,'Fontsize',20)

% Relative density calculation
Relative_matrix = ((VolumeObject - sum(Data_matrix(:,6)))/VolumeObject)*100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Prediction part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% confusion result and model accuracy
%load the prdiction model result
load('X:\Castro\recon_pore\matfile\July15result.mat'); % < - dim_3D is mat file

%% data cleaning (duplicate data has to be eliminated)
dim_3D = sortrows(unique(dim_3D(:,:),'rows'),3);

% change the scale of the datapoint data reconstruction
%% data cleaning (based on the coorelation cofficient {0.8~0.95} in order to preserve the data and eliminate )
dim_3D(dim_3D(:,4)>0.95,:) = [];
dim_3D(dim_3D(:,4)< 0.8,:) = [];

% sort the data by height
Data_matrix_mesh = sortrows(Data_matrix_mesh,3);

%% Confusion matrix and Final result
%%%%%%%%%%%%%%%%%%%%%%%%%%%Confusion Matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the idea of confusion matrix is padding and the confusion matrix has to be improved

range = 10; %The Padd_path
% 
TP = 0; % True positive in any direction
TNX = 0;% True negative in x direction
TNY = 0;% True negative in y direction
TNZ = 0;% True negative in z direction
% 
for index_CT = 1: size(Data_matrix_mesh,1)
     for height = 1:max(dim_3D(:,3))
         for row = 1: max(dim_3D(:,2)) % y_dir
             for col = 1:max(dim_3D(:,1)) % x_dir
                 % searching the the data point around twenty point (large ragne searching)
                 if ~(height <= range || row <= range || col <=range) && ~(height >(max(dim_3D(:,3)) - range )|| row >(max(dim_3D(:,2))-range) || (col > max(dim_3D(:,1))-range) )
                     if Data_matrix_mesh(index_CT,3) == height
                         for index_pre = 1: size(dim_3D,1)
                            if dim_3D(index_pre,3) <= (height + range) && dim_3D(index_pre,3) >= (height - range)
                                if dim_3D(index_pre,2) <= (row + range) && dim_3D(index_pre,2) >= (row - range)
                                    if dim_3D(index_pre,1) <= (col+range) && dim_3D(index_pre,1) >= (col - range)
                                        TP = TP+1; %TP
                                    else
                                        TNX = TNX+1; % TN
                                    end
                                else
                                    TNY = TNY + 1;%TN
                                end
                            else
                                TNZ = TNZ+1; %TN
                            end
                         end
                     end
                 end
             end
         end
     end
 end
 
%% due to the result has repeated so that it has to be divide  by the cycle number
% The cycle number is size of Data_matrix_mesh
TP = round(TP/size(Data_matrix_mesh,1));
TN = round((TNX+TNY+TNY)/size(Data_matrix_mesh,1));
FN = range*range*range-(range*3);
FP = range*range*range;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the new pore volume and location
pre_matrix = zeros(max(dim_3D(:,1)),max(dim_3D(:,2)),max(dim_3D(:,3)));
dim_3D_test = dim_3D;
%%%%%%%%%%%%%%%%%%%% Using padding matrix re-arange the pore data location %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Still testing  %%%%The new_3D matrix has to be exsitence%%%%%%%%%%%

% kernel size
k_x = 31; % x-dir
k_y = 31; % y-dir
k_z = 27; % z-dir

%% "new"  calculate the new pore volume and location

for index_dim3D = 1:size(dim_3D_test,1) % searching index of pore data
    for height = 1:max(dim_3D_test(:,3))
        if height == dim_3D_test(index_dim3D,3)
            for row = 1:max(dim_3D_test(:,2))
                if row == dim_3D_test(index_dim3D,2)
                    for col = 1:max(dim_3D_test(:,1))
                        if col == dim_3D_test(index_dim3D,1)
                            pre_matrix(col,row,height) = dim_3D_test(index_dim3D,4);
                        else
                            continue
                        end
                    end
                else
                	continue
                end
            end
        else
            continue
        end
    end
end

% make a table as same as dim_3D x,y,z,value
new_3D = [];
index = 1;
for height = 1: size(pre_matrix,3)
    for row = 1:size(pre_matrix,2)
        for col = 1:size(pre_matrix,1)
            if pre_matrix(col,row,height) ~= 0
                new_3D(index,1) = col;
                new_3D(index,2) = row;
                new_3D(index,3) = height;
                new_3D(index,4) = pre_matrix(col,row,height);
                index = index +1;
            end
        end
    end
end

if isempty(new_3D) == 1
    new_3D(1,1) = 0;
    new_3D(1,2) = 0;
    new_3D(1,3) = 0;
    new_3D(1,4) = 1;
end
    
figure, scatter3(new_3D(:,1),new_3D(:,2),new_3D(:,3),40,new_3D(:,4),'filled')
title('pre_matrix_first')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
zlim([0,120])
set(gca,'Fontsize',20)
%% Padding and Creating new pore structure based upon Post-processing
% kernel size (have to chose the odd number)
k_x = 31;
k_y = 21;
k_z = 3;

for height = 1: size(pre_matrix,3)
    for row = 1: size(pre_matrix,2)
        for col = 1:size(pre_matrix,1)
            if (height > ((k_z-1)/2)) && (height < size(pre_matrix,3)-((k_z-1)/2)) && (row > ((k_y-1)/2)) && (row < size(pre_matrix,2)-((k_y-1)/2)) && (col > ((k_x-1)/2)) && (col < size(pre_matrix,1)-((k_x-1)/2))
                %while pre_matrix(col,row,height)~= 0
                    for kernel_z = -((k_z-1)/2):((k_z-1)/2)
                        for kernel_y = -((k_y-1)/2):((k_y-1)/2)
                            for kernel_x = -((k_x-1)/2):((k_x-1)/2)
                                if ((height+kernel_z)~=height) && ((row+kernel_y) ~= row) && ((col+kernel_x) ~= col)
                                    pre_matrix(col,row,height) = pre_matrix(col,row,height) + pre_matrix(col+kernel_x,row+kernel_y,height+kernel_z);
                                    pre_matrix(col+kernel_x,row+kernel_y,height+kernel_z) = 0;
                                end
                            end
                        end
                    end 
                %end
            end
        end
    end
end

new_3D = [];
index = 1;
for height = 1: size(pre_matrix,3)
    for row = 1:size(pre_matrix,2)
        for col = 1:size(pre_matrix,1)
            if pre_matrix(col,row,height) ~= 0
                new_3D(index,1) = col;
                new_3D(index,2) = row;
                new_3D(index,3) = height;
                new_3D(index,4) = pre_matrix(col,row,height);
                index = index +1;
            end
        end
    end
end

if isempty(new_3D) == 1
    new_3D(1,1) = 0;
    new_3D(1,2) = 0;
    new_3D(1,3) = 0;
    new_3D(1,4) = 1;
end
    
figure, scatter3(new_3D(:,1),new_3D(:,2),new_3D(:,3),40,new_3D(:,4),'filled')
title('new_3D')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
zlim([0,120])
set(gca,'Fontsize',20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the original prediction result in 3D 
figure, scatter3(dim_3D(:,1),dim_3D(:,2),dim_3D(:,3),40,dim_3D(:,4),'filled')
title('dim_3D')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
zlim([0,120])
set(gca,'Fontsize',20)

%% the pore size of prediction and the pore size of CT (Ratio)
%% The Tunning part %%
new_3D_scale = [];
new_3D_scale = new_3D(:,4)*80;


%% Plot the X-ray and prediction model 
% change the x,y,z label to the same scale (x,y) in micrometer z in layer
Voxel_size = 15.245;
%% prediction model

figure, scatter3(new_3D(:,1)*Voxel_size,new_3D(:,2)*Voxel_size,new_3D(:,3),40,new_3D_scale,'filled')
title('Prediction')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
zlim([0,120])
set(gca,'Fontsize',20)

%% Original Prediction model
figure, scatter3(dim_3D(:,1)*Voxel_size,dim_3D(:,2)*Voxel_size,dim_3D(:,3),40,dim_3D(:,4),'filled')
title('Original Prediction')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
% zlim([0,100])
set(gca,'Fontsize',20)
%% CT_result
figure, scatter3(Data_matrix_mesh(:,1)*Voxel_size,Data_matrix_mesh(:,2)*Voxel_size,Data_matrix_mesh(:,3),40,Data_matrix_mesh(:,4),'filled')
title('CT model')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
% zlim([0,100])
% xlim([0, 4500])
% ylim([0, 2500])
set(gca,'Fontsize',20)

%%  The perfoemace of model
Accuracy = (TP/(TP+FP))*100; %% Precision of model
Recall = (TP/(TP+FN))*100;  %% sensitive of the model
F_score = (2*(Accuracy*Recall))/(Accuracy+Recall); %% the mono indication => entirely indication!!

%% Denstiy calculation:
% Object dimension (micro meter)
H = 2000; 
L = 4000;
W = 2000;
Object_volume = H*L*W;
% CT Porosity
CT_density = ((Object_volume-sum(Data_matrix_mesh(:,4)))/Object_volume)*100;
Pre_density = ((Object_volume-sum(new_3D_scale))/Object_volume)*100;
%%  display the result on the matlab command window
disp(['_______The result of comparison________']);
disp(['TP:' num2str(TP) ])
disp(['TN:' num2str(TN) ])
disp(['FP:' num2str(FP) ])
disp(['FN:' num2str(FN) ])
disp(['Accuracy:' num2str(Accuracy) '%'])
disp(['Recall:' num2str(Recall) '% (Sensitive)'])
disp(['F_score:' num2str(F_score) '%']);
disp (['Object Volume:' num2str(Object_volume) 'um^3'])
disp(['CT Result_______________________'])
disp(['CT_pore_volume:' num2str(sum(Data_matrix_mesh(:,4))) 'um^3'])
disp(['CT density:' num2str(CT_density) '%'])
disp(['Prediction Resul________________'])
disp(['Pre_pore_volume:' num2str(sum(new_3D_scale)) 'um^3'])
disp(['Pre_densityt' num2str(Pre_density) '%'])
