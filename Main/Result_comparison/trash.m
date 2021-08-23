% % clc,clear
% % close all
% % load('prediction_result.mat')
% % load('ct_result.mat')
% % dim_3D = sortrows(unique(dim_3D(:,:),'rows'),3);
% % 
% % dim_3D_ch = dim_3D; % duplicate the matrix
% % dim_3D_new = [];
% % nn = 1; %
% % 
% % 
% % % figure, scatter3(dim_3D(:,1),dim_3D(:,2),dim_3D(:,3),40,dim_3D(:,4),'filled')
% % % xlabel('X-axis(\mum)')
% % % ylabel('Y-axis(\mum)')
% % % zlabel('Layer')
% % % colorbar('East')
% % % zlim([0,120])
% % % set(gca,'Fontsize',20)
% % 
% % 
% % 
% % 
% % % read the data from dim_3D(prediction model)
% % % build the 3-D matrix
% % pre_matrix = zeros(max(dim_3D(:,1)),max(dim_3D(:,2)),max(dim_3D(:,3)));
% % for z_dir = 1: size(pre_matrix,3)
% %     for y_dir = 1: size(pre_matrix,2)
% %         for x_dir = 1: size(pre_matrix,1)
% %             for ii = 1: size(dim_3D,1)
% %                 if ((x_dir == dim_3D(ii,1)) && (y_dir == dim_3D(ii,2)) && (z_dir == dim_3D(ii,3)))
% %                     pre_matrix(x_dir,y_dir,z_dir) = dim_3D(ii,4);
% %                 end
% %             end
% %         end
% %     end
% % end
% % 
% % % read the data from Data_matrix_mesh (CT result)
% % CT_model = zeros(max(Data_matrix_mesh(:,1)),max(Data_matrix_mesh(:,2)),max(Data_matrix_mesh(:,3)));
% % for z_dir = 1: size(CT_model,3)
% %     for y_dir = 1: size(CT_model,2)
% %         for x_dir = 1: size(CT_model,1)
% %             for ii = 1: size(Data_matrix_mesh,1)
% %                 if ((x_dir == Data_matrix_mesh(ii,1)) && (y_dir == Data_matrix_mesh(ii,2)) && (z_dir == Data_matrix_mesh(ii,3)))
% %                     CT_model(x_dir,y_dir,z_dir) = Data_matrix_mesh(ii,4);
% %                 end
% %             end
% %         end
% %     end
% % end
% % 
% % 
% % 
% % %............................................................................................
% % % objec-> kernel -> moving
% % 
% % % kernel size
% k_x = 31; % x-dir
% k_y = 31; % y-dir
% k_z = 27; % z-dir
% 
% % for height = 1:(max(dim_3D(:,3)))
% %     for row = 1:(max(dim_3D(:,2)))
% %         for col = 1:(max(dim_3D(:,1)))
% %             % build the kernel 
% %             % boundary problem
% %             if (height > ((k_z-1)/2)) && (height < max(dim_3D(:,3))-((k_z-1)/2)) && (row > ((k_y-1)/2)) && (row < max(dim_3D(:,2))-((k_y-1)/2)) && (col > ((k_x-1)/2)) && (col < max(dim_3D(:,1))-((k_x-1)/2))
% %                 % where is the original pore location, it is a center of
% %                 % kenrnal
% %                 if pre_matrix(col,row,height) ~= 0
% %                    for k = -((k_z-1)/2):((k_z-1)/2)
% %                        for j = -((k_y-1)/2):((k_y-1)/2)
% %                            for i = -((k_x-1)/2):((k_x-1)/2)
% %                                % add the operation
% %                                if ((height+k)==height) && ((row+j) == row) && ((col+i) == col)
% %                                   continue
% %                                else
% %                                    disp(['coorinnate;' num2str(col) ' ' num2str(row) ' ' num2str(height)])
% %                                    disp(['original' num2str( pre_matrix(col,row,height))])
% %                                    disp(['next' num2str(pre_matrix(col+i,row+j,height+k))])
% %                                    pre_matrix(col,row,height) = pre_matrix(col+i,row+j,height+k)+pre_matrix(col,row,height);
% %                                    pre_matrix(col+i,row+j,height+k) = 0; 
% %                                    disp(['after' num2str(pre_matrix(col,row,height))])
% %                                    disp(['...............'])                    
% %                                 end
% %                             end
% %                         end
% %                     end
% %                 end
% %             else
% %                 % Due to this region was neglect there will be zero and
% %                 % this is boundary
% %                 pre_matrix(col,row,height) = 0;
% %             end
% %         end
% %     end
% % end
% % % 
% % % 
% % % %% buit the new confusion matrix in order to calculate the real number and data visualization code
% % % TP = 0;
% % % FN = 0;
% % % FP = 0;
% % % TN = 0;
% % % % compare two matrix 1.pre_matrix 2, CT_model
% % % k_x = 9;
% % % k_y = 9;
% % % k_z = 3;
% % % for height = 1:(max((Data_matrix_mesh(:,3))))
% % %     for row = 1:(max(dim_3D(:,2)))
% % %         for col = 1:(max(Data_matrix_mesh(:,1)))
% % %             % set the kernel and searching
% % %             % avoid the bounday problem 
% % %             if (height > ((k_z-1)/2)) && (height < max((Data_matrix_mesh(:,3)))-((k_z-1)/2)) && (row > ((k_y-1)/2)) && (row < max(dim_3D(:,2))-((k_y-1)/2)) && (col > ((k_x-1)/2)) && (col < max((Data_matrix_mesh(:,3)))-((k_x-1)/2))
% % %                 % the main zone
% % %                 if CT_model(col,row,height) ~= 0   % actula yes
% % %                      for k = -((k_z-1)/2):((k_z-1)/2)
% % %                        for j = -((k_y-1)/2):((k_y-1)/2)
% % %                            for i = -((k_x-1)/2):((k_x-1)/2)
% % %                                % searching the cubic based on the keneal 
% % %                                if pre_matrix(col+i,row+j,height+k) ~=0 % predict yes
% % %                                    disp(['the prediction_point:' num2str(col+i) ' ' num2str(row+j) ' ' num2str(height+k)])
% % %                                    disp(['the CT_point:' num2str(col) ' ' num2str(row) ' ' num2str(height)])
% % %                                    TP = TP+1;
% % %                                else % predict NO
% % %                                    FN = FN+1;
% % %                                end
% % %                             end
% % %                         end
% % %                      end
% % %                      FN = FN - (i*j*k);
% % %                 elseif CT_model(col,row,height) == 0 % Actual No
% % %                     for k = -((k_z-1)/2):((k_z-1)/2)
% % %                        for j = -((k_y-1)/2):((k_y-1)/2)
% % %                            for i = -((k_x-1)/2):((k_x-1)/2)
% % %                                % searching the cubic based on the keneal 
% % %                                if pre_matrix(col+i,row+j,height+k) ~=0 % predict yes
% % %                                    FP = FP+1;
% % %                                else % predict NO
% % %                                    TN = TN+1;
% % %                                end
% % %                             end
% % %                         end
% % %                     end
% % %                     FP = FP - (i*k*j);
% % %                     TN = TN - (i*k*j);
% % %                 else
% % %                     continue;
% % %                 end
% % %             end
% % %         end
% % %     end
% % % end
% % % % the total point:
% % % total_point = max(Data_matrix_mesh(:,1))*max(dim_3D(:,2))*max(Data_matrix_mesh(:,3));
% % % disp(['Totoal point:' num2str(total_point)])
% % % disp(['TP:' num2str(TP) ])
% % % disp(['TN:' num2str(TN) ])
% % % disp(['FP:' num2str(FP) ])
% % % disp(['FN:' num2str(FN) ])
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % make a table as same as dim_3D x,y,z,value
% % % new_3D = [];
% % % index = 1;
% % % for height = 1: size(pre_matrix,3)
% % %     for row = 1:size(pre_matrix,2)
% % %         for col = 1:size(pre_matrix,1)
% % %             if pre_matrix(col,row,height) ~= 0
% % %                 new_3D(index,1) = col;
% % %                 new_3D(index,2) = row;
% % %                 new_3D(index,3) = height;
% % %                 new_3D(index,4) = pre_matrix(col,row,height);
% % %                 index = index +1;
% % %             end
% % %         end
% % %     end
% % % end
% % % % 
% % % figure, scatter3(new_3D(:,1),new_3D(:,2),new_3D(:,3),40,new_3D(:,4),'filled')
% % % xlabel('X-axis(\mum)')
% % % ylabel('Y-axis(\mum)')
% % % zlabel('Layer')
% % % colorbar('East')
% % % zlim([0,120])
% % % set(gca,'Fontsize',20)
% % % 
% % % figure, scatter3(Data_matrix_mesh(:,1),Data_matrix_mesh(:,2),Data_matrix_mesh(:,3),40,Data_matrix_mesh(:,4),'filled')
% % % xlabel('X-axis(\mum)')
% % % ylabel('Y-axis(\mum)')
% % % zlabel('Layer')
% % % colorbar('East')
% % % zlim([0,120])
% % % set(gca,'Fontsize',20)
% % 
% % %% debug code
% % a = max(Data_matrix_mesh(:,4));
% % b = min(Data_matrix_mesh(:,4));
% % c = max(new_3D(:,4));
% % d = min(new_3D(:,4));
% % % In the same scale
% % Nominal_new_3D = (new_3D(:,4)-d)/(c-d);
% % % new_3D_scale = (Nominal_new_3D * (a-b))+b;
% % new_3D_scale = Nominal_new_3D*160;
% % %% Plot the X-ray and prediction model 
% % % change the x,y,z label to the same scale (x,y) in micrometer z in layer
% % Voexl_size = 18.605;
% % %% prediction model
% % 
% % figure, scatter3(new_3D(:,1)*Voxel_size,new_3D(:,2)*Voxel_size,new_3D(:,3),40,new_3D_scale,'filled')
% % xlabel('X-axis(\mum)')
% % ylabel('Y-axis(\mum)')
% % zlabel('Layer')
% % colorbar('East')
% % zlim([0,120])
% % set(gca,'Fontsize',20)
% % 
% % %% Original Prediction model
% % figure, scatter3(dim_3D(:,1)*Voxel_size,dim_3D(:,2)*Voxel_size,dim_3D(:,3),40,dim_3D(:,4),'filled')
% % title('Original Prediction')
% % xlabel('X-axis(\mum)')
% % ylabel('Y-axis(\mum)')
% % zlabel('Layer')
% % colorbar('East')
% % zlim([0,120])
% % set(gca,'Fontsize',20)
% % %% CT_result
% % figure, scatter3(Data_matrix_mesh(:,1)*Voxel_size,Data_matrix_mesh(:,2)*Voxel_size,Data_matrix_mesh(:,3),40,Data_matrix_mesh(:,4),'filled')
% % xlabel('X-axis(\mum)')
% % ylabel('Y-axis(\mum)')
% % zlabel('Layer')
% % colorbar('East')
% % zlim([0,120])
% % set(gca,'Fontsize',20)
% % 
% % %%  The perfoemace of model
% % Accuracy = (TP/(TP+FP))*100; %% Precision of model
% % Recall = (TP/(TP+FN))*100;  %% sensitive of the model
% % F_score = (2*(Accuracy*Recall))/(Accuracy+Recall); %% the mono indication => entirely indication!!
% % 
% % %% Denstiy calculation:
% % % Object dimension (micro meter)
% % H = 2000; 
% % L = 4000;
% % W = 2000;
% % Object_volume = H*L*W;
% 
% % % CT Porosity
% CT_density = ((Object_volume-sum(Data_matrix_mesh(:,4)))/Object_volume)*100;
% Pre_density = ((Object_volume-sum(new_3D_scale))/Object_volume)*100;
% 
% %%  display the result on the matlab command window
% % disp(['_______The result of comparison________']);
% % disp(['TP:' num2str(TP) ])
% % disp(['TN:' num2str(TN) ])
% % disp(['FP:' num2str(FP) ])
% % disp(['FN:' num2str(FN) ])
% % disp(['Accuracy:' num2str(Accuracy) '%'])
% % disp(['Recall:' num2str(Recall) '% (Sensitive)'])
% % disp(['F_score:' num2str(F_score) '%']);
% % disp (['Object Volume:' num2str(Object_volume) 'um^3'])
% % disp(['CT Result_______________________'])
% % disp(['CT_pore_volume:' num2str(sum(Data_matrix_mesh(:,4))) 'um^3'])
% % disp(['CT density:' num2str(CT_density) '%'])
% % disp(['Prediction Resul________________'])
% % disp(['Pre_pore_volume:' num2str(sum(new_3D_scale)) 'um^3'])
% % disp(['Pre_densityt' num2str(Pre_density) '%'])
% 
%%  the code desin( Confusion matrix, The matrix)
figure, scatter3(dim_3D(1:10,1)*Voxel_size,dim_3D(1:10,2)*Voxel_size,dim_3D(1:10,3),40,dim_3D(1:10,4),'filled')    
figure, scatter3(Data_matrix_mesh(1:10,1),Data_matrix_mesh(1:10,2),Data_matrix_mesh(1:10,3),40,Data_matrix_mesh(1:10,4),'filled')























