clc,clear
close all
% Test 82
% Test 83
path = 'X:\Castro\pore_prediction_NOS\Test21.xlsx';
[~,sheet_name]=xlsfinfo(path);
for k=1:numel(sheet_name)
  data_sheet{k}=xlsread(path,sheet_name{k});
end

Meltpool_dimension = data_sheet{1,1}(1:size(data_sheet{1,1},1),1:size(data_sheet{1,1},2))'; %16x 1 Table
%Number_of_Meltpool = data{1,1}(1:16,1:3); % 16 x1 Table

corr_melt = corr(Meltpool_dimension);
figure,surf(corr_melt),colorbar   %/max(max(a(1:16,1:16)))
for i = 1:size(corr_melt,1)
    med(i) = median(corr_melt(:,i));
end