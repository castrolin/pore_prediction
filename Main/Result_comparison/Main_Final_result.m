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
%% read the file
% x = xlsread('Pore.xlsx');
x = xlsread('X:\Castro\recon_pore\excel\Pore_fist_July15.xlsx');
%% Parameter setup
Voxel_size = 18.604;
meshsize = 0.018604; % 20um 0.02 mm
% object_dimension (micrometer)
Height = 2000;
Width  = 3000;
Long = 4000;
VolumeObject = Height * Width * Long;
%% creat the data matrix
diameter = x(:,3);
voxel = x(:,8);
volume = x(:,7);

%% the coordinate  (different sample has differnet cooridnate !!)
y_dir = x(:,6); % soft -> software
z_dir = x(:,4);
x_dir = x(:,5);

% Data_matrix = [round(x_dir/meshsize),round(y_dir/meshsize),round(z_dir/meshsize),voxel,diameter];z->x, y->z, x->y
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


%% confusion result and model accuracy

%load the prdiction model result

% load('prediction_result.mat');
% load('pre_spatter_zero.mat');
% load('pre_spatter_0to30.mat');
% load('pre_spatter_30to60.mat');
% load('pre_spatter_60to90.mat');
% load('thesis_20deg.mat');  % < - dim_3D is mat file
load('X:\Castro\recon_pore\matfile\July15result.mat');

%% data cleaning (duplicate data has to be eliminated)

dim_3D = sortrows(unique(dim_3D(:,:),'rows'),3);

%% dimenstionto mesh 
Data_matrix_mesh = sortrows([round(Data_matrix(:,1)/Voxel_size),round(Data_matrix(:,2)/Voxel_size),round(Data_matrix(:,3)/35)+1,Data_matrix(:,5)],3);
% change the scale of the datapoint data reconstruction
%% data cleaning
%Data_matrix_mesh(Data_matrix_mesh(:,1)<50,:) = []; % based on the x-dirction fit the scale of prediction 107
%Data_matrix_mesh(Data_matrix_mesh(:,1)>300,:) = []; % based on the y-dirction fit the scale of prediction 296 (coordinate)

dim_3D(dim_3D(:,4)>0.95,:) = [];
dim_3D(dim_3D(:,4)< 0.8,:) = [];
%% rescaled
% OriMin = min(Data_matrix_mesh(:,1)); % original min
% OriMax = max(Data_matrix_mesh(:,1)); % original max
% Data_matrix_mesh(:,1) = Data_matrix_mesh(:,1) - 50 +1; %- 107 + 1;
% sort the data by height
Data_matrix_mesh = sortrows(Data_matrix_mesh,3);


%% Confusion matrix and Final result
% the scale i s wrong

range = 10; %The Padd_path
% 
TP = 0;
TNX = 0;
TNY = 0;
TNZ = 0;
% 
% for index_CT = 1: size(Data_matrix_mesh,1)
%     for height = 1:max(dim_3D(:,3))
%         for row = 1: max(dim_3D(:,2)) % y_dir
%             for col = 1:max(dim_3D(:,1)) % x_dir
%                 % searching the the data point around twenty point (large ragne searching)
%                 if ~(height <= range || row <= range || col <=range) && ~(height >(max(dim_3D(:,3)) - range )|| row >(max(dim_3D(:,2))-range) || (col > max(dim_3D(:,1))-range) )
%                     if Data_matrix_mesh(index_CT,3) == height
%                         for index_pre = 1: size(dim_3D,1)
%                            if dim_3D(index_pre,3) <= (height + range) && dim_3D(index_pre,3) >= (height - range)
%                                if dim_3D(index_pre,2) <= (row + range) && dim_3D(index_pre,2) >= (row - range)
%                                    if dim_3D(index_pre,1) <= (col+range) && dim_3D(index_pre,1) >= (col - range)
%                                        TP = TP+1; %TP
%                                    else
%                                        TNX = TNX+1; % TN
%                                    end
%                                else
%                                    TNY = TNY + 1;%TN
%                                end
%                            else
%                                TNZ = TNZ+1; %TN
%                            end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

%% calculate the new pore volume and location
pre_matrix = zeros(max(dim_3D(:,1)),max(dim_3D(:,2)),max(dim_3D(:,3)));
dim_3D_test = dim_3D;

% kernel size
k_x = 31; % x-dir
k_y = 31; % y-dir
k_z = 27; % z-dir


% objec-> kernel -> moving
% for height = 1:(max(dim_3D(:,3)))
%     for row = 1:(max(dim_3D(:,2)))
%         for col = 1:(max(dim_3D(:,1)))
%             % build the kernel 
%             % boundary problem (avoid boundary)
%             if (height > ((k_z-1)/2)) && (height < max(dim_3D(:,3))-((k_z-1)/2)) && (row > ((k_y-1)/2)) && (row < max(dim_3D(:,2))-((k_y-1)/2)) && (col > ((k_x-1)/2)) && (col < max(dim_3D(:,1))-((k_x-1)/2))
%                 % where is the original pore location, it is a center of
%                 % kenrnal
%                 if pre_matrix(col,row,height) ~= 0
%                    for k = -((k_z-1)/2):((k_z-1)/2)
%                        for j = -((k_y-1)/2):((k_y-1)/2)
%                            for i = -((k_x-1)/2):((k_x-1)/2)
%                                % add the operation
%                                if ((height+k)==height) && ((row+j) == row) && ((col+i) == col)
%                                   continue
%                                else
%                                    disp(['coorinnate;' num2str(col) ' ' num2str(row) ' ' num2str(height)])
%                                    disp(['original' num2str( pre_matrix(col,row,height))])
%                                    disp(['next' num2str(pre_matrix(col+i,row+j,height+k))])
%                                    pre_matrix(col,row,height) = pre_matrix(col+i,row+j,height+k)+pre_matrix(col,row,height);
%                                    pre_matrix(col+i,row+j,height+k) = 0; % make the be sum point become a zero
%                                    disp(['after' num2str(pre_matrix(col,row,height))])
%                                    disp(['...............'])                    
%                                 end
%                             end
%                         end
%                     end
%                 end
%             else
%                 % Due to this region was neglect there will be zero and
%                 % this is boundary
%                 pre_matrix(col,row,height) = 0;
%             end
%         end
%     end
% end


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
% make a table as same as dim_3D x,y,z,value
% new_3D = [];
% index = 1;
% for height = 1: size(pre_matrix,3)
%     for row = 1:size(pre_matrix,2)
%         for col = 1:size(pre_matrix,1)
%             if pre_matrix(col,row,height) ~= 0
%                 new_3D(index,1) = col;
%                 new_3D(index,2) = row;
%                 new_3D(index,3) = height;
%                 new_3D(index,4) = pre_matrix(col,row,height);
%                 index = index +1;
%             end
%         end
%     end
% end

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
% make a table as same as dim_3D x,y,z,value
% new_3D = [];
% index = 1;
% for height = 1: size(pre_matrix,3)
%     for row = 1:size(pre_matrix,2)
%         for col = 1:size(pre_matrix,1)
%             if pre_matrix(col,row,height) ~= 0
%                 new_3D(index,1) = col;
%                 new_3D(index,2) = row;
%                 new_3D(index,3) = height;
%                 new_3D(index,4) = pre_matrix(col,row,height);
%                 index = index +1;
%             end
%         end
%     end
% end

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

figure, scatter3(dim_3D(:,1),dim_3D(:,2),dim_3D(:,3),40,dim_3D(:,4),'filled')
title('dim_3D')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
zlim([0,120])
set(gca,'Fontsize',20)
%% due to the result has repeated so that it has to be divide  by the cycle number
% The cycle number is size of Data_matrix_mesh
TP = round(TP/size(Data_matrix_mesh,1));
TN = round((TNX+TNY+TNY)/size(Data_matrix_mesh,1));
FN = range*range*range-(range*3);
FP = range*range*range;
%% the pore size of prediction and the pore size of CT (Ratio)
% in order to have same scale bar!
a = max(Data_matrix_mesh(:,4));
b = min(Data_matrix_mesh(:,4));
c = max(new_3D(:,4));
d = min(new_3D(:,4));
% In the same scale
Nominal_new_3D = (new_3D(:,4)-d)/(c-d);
new_3D_scale = [];
% for lengthofNominal_new_3D = 1:size(Nominal_new_3D,1)
%     new_3D_scale(lengthofNominal_new_3D,1) = Nominal_new_3D(lengthofNominal_new_3D,1) *(a-b)+b;
% end
new_3D_scale = new_3D(:,4)*80;
% new_3D_scale = (Nominal_new_3D .* (a-b))+b;
% new_3D_scale = Nominal_new_3D;
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