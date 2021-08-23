%% confusion result and model accuracy

%load the 
load('prediction_result.mat');

%% data cleaning (duplicate data has to be eliminated)

dim_3D = sortrows(unique(dim_3D(:,:),'rows'),3);

%% dimenstionto mesh 
Data_matrix_mesh = sortrows([round(Data_matrix(:,1)/18.605),round(Data_matrix(:,2)/18.605),round(Data_matrix(:,3)/18.605)+1,Data_matrix(:,5)],3);
% change the scale of the datapoint data reconstruction
% data cleaning
Data_matrix_mesh(Data_matrix_mesh(:,1)<107,:) = [];
Data_matrix_mesh(Data_matrix_mesh(:,1)>296,:) = [];
% rescaled
OriMin = min(Data_matrix_mesh(:,1)); % original min
OriMax = max(Data_matrix_mesh(:,1)); % original max
Data_matrix_mesh(:,1) = Data_matrix_mesh(:,1) - 107 + 1;
% sort the data by height
Data_matrix_mesh = sortrows(Data_matrix_mesh,3);

%% the pore size of prediction and the pore size of CT (Ratio)
a = max(Data_matrix_mesh(:,4));
b = min(Data_matrix_mesh(:,4));
c = max(dim_3D(:,4));
d = min(dim_3D(:,4));
% In the same scale
Nominal_dim_3D = (dim_3D(:,4)-d)/(c-d);
dim_3D_scale = (Nominal_dim_3D * (a-b))+b;

%% Confusion matrix and Final result
% the scale i s wrong
range = 2; %The Padd_path
TP = 0;
TNX = 0;
TNY = 0;
TNZ = 0;

for index_CT = 1: size(Data_matrix_mesh,1)
    for height = 1:max(dim_3D(:,3))
        for row = 1: max(dim_3D(:,2)) % y_dir
            for col = 1:max(dim_3D(:,1)) % x_dir
                % searching the the data point around twenty point (large ragne searching)
                if ~(height <= range || row <= range || col <=range) || ~(height >(max(dim_3D(:,3)) - range )|| row >(max(dim_3D(:,2))-range) || (col > max(dim_3D(:,1))-range) )
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
%% Plot the X-ray and prediction model 
% change the x,y,z label to the same scale (x,y) in micrometer z in layer
Voexl_size = 18.605;
%% prediction model

figure, scatter3(dim_3D(:,1)*Voxel_size,dim_3D(:,2)*Voxel_size,dim_3D(:,3),40,dim_3D_scale,'filled')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
zlim([0,120])

%% CT_result
figure, scatter3(Data_matrix_mesh(:,1)*Voxel_size,Data_matrix_mesh(:,2)*Voxel_size,Data_matrix_mesh(:,3),40,Data_matrix_mesh(:,4),'filled')
xlabel('X-axis(\mum)')
ylabel('Y-axis(\mum)')
zlabel('Layer')
colorbar('East')
zlim([0,120])


%%  The perfoemace of model
Accuracy = (TP/(TP+FP))*100; %% Precision of model
Recall = (TP/(TP+FN))*100;  %% sensitive of the model
F_score = (2*(Accuracy*Recall))/(Accuracy+Recall); %% the mono indication => entirely indication!!

%% Denstiy calculation:
% Object dimension (micro meter)
H = 2000; 
L = 3500;
W = 2000;
Object_volume = H*L*W;
% CT Porosity
CT_density = ((Object_volume-sum(Data_matrix_mesh(:,4)))/Object_volume)*100;
Pre_density = ((Object_volume-sum(dim_3D_scale))/Object_volume)*100;
%%  display the result on the matlab command window
disp(['_______The result of comparison________']);
disp(['Accuracy:' num2str(Accuracy) '%'])
disp(['Recall:' num2str(Recall) '% (Sensitive)'])
disp(['F_score:' num2str(F_score) '%']);
disp (['Object Volume:' num2str(Object_volume) 'um^3'])
disp(['CT Result_______________________'])
disp(['CT_pore_volume:' num2str(sum(Data_matrix_mesh(:,4))) 'um^3'])
disp(['CT density:' num2str(CT_density) '%'])
disp(['Prediction Resul________________'])
disp(['Pre_pore_volume:' num2str(sum(dim_3D_scale)) 'um^3'])
disp(['Pre_densityt' num2str(Pre_density) '%'])