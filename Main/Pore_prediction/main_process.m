%%%%% Full structure of pore prediction structure
%%% Combine the Analysis and prediction method
%%% This is the final code
%% read viedo or file
%% SOM algorithm
clc,clear
close all

% files = {'Normal1.mat','Normal2.mat','Normal3.mat',...
%         'Normal4.mat','Normal5.mat','Normal6.mat','Normal7.mat',...
%         'Abnormal1.mat','Abnormal2.mat','Abnormal3.mat',...
%         'Abnormal4.mat','Abnormal5.mat','Abnormal6.mat'};

n = 1 ;
name = {};
for datanum = 1:100 % video name of file
    files{n} = strcat(sprintf('%1.f',datanum),'.mat');
    n = n+1;
end
%% read the Numner of cluster excel in order to define K-value
NOC_excel = xlsread("Number_of_cluster_0715.xlsx"); % NOC = Number Of Cluster

% centromatrix size => Number of Output X Number of Input(width length
% ratio) 3*3 ouput
% 7X7
%% testing the optimize cluster in code
%data = load(files{1});
 % record the all of the centorid(Number of data,Number of cluster,Number of feature)
Layer = 1;
for i = 1:numel(files)
    worksheetName = sprintf('Test%d.xlsx',i) ;
       %% post-pocessing 
    % for each files load and setup
    Data = load(files{i});
    data_set = [[Data.DataBase.Length];[Data.DataBase.Width];[Data.DataBase.Ratio];[Data.DataBase.Angle];[Data.DataBase.NOS]]';
    
    %[IDX,C,SUMD,K]=best_kmeans(data_set);
    K = NOC_excel(i,2);
    % if cluseter number is equal to one (the correlation will failed)
    if K == 1
        K = 2;
    end
    centromatrix = zeros(numel(files),K,5);
    %-------------------------------------
    leng(1,1:size([Data.DataBase(1:end)],2))= [Data.DataBase(1:end).Length];
    width(1,1:size([Data.DataBase(1:end)],2))= [Data.DataBase(1:end).Width];
    ratio(1,1:size([Data.DataBase(1:end)],2))= [Data.DataBase(1:end).Ratio];
    angle(1,1:size([Data.DataBase(1:end)],2))= [Data.DataBase(1:end).Angle];
    nos(1,1:size([Data.DataBase(1:end)],2))= [Data.DataBase(1:end).NOS];
   n=1;
   for k = 1:size(leng,1)
       for j = 1:size(leng,2)
           if ratio(k,j)~=0
               lengthdata(n)=leng(k,j);
               n= n+1;
           end
       end
   end
   
   n=1;
   for k = 1:size(width,1)
       for j = 1:size(width,2)
           if ratio(k,j)~=0
               widthdata(n)=width(k,j);
               n= n+1;
           end
       end
   end
   
   n=1;
   for k = 1:size(ratio,1)
       for j = 1:size(ratio,2)
           if ratio(k,j)~=0
               ratiodata(n)=ratio(k,j);
               n= n+1;
           end
       end
   end
   n =1;
   for k = 1:size(angle,1)
       for j = 1:size(angle,2)
           if ratio(k,j)~=0   % in order to keep same size as other array
               angledata(n)=angle(k,j);
               n= n+1;
           end
       end
   end
   
   n =1;
   for k = 1:size(nos,1)
       for j = 1:size(nos,2)
           if ratio(k,j)~=0   % in order to keep same size as other array
               nosdata(n)=nos(k,j);
               n= n+1;
           end
       end
   end

   %% Self-organize map running
   if exist('lengthdata','var') ~= 0 || exist('widthdata','var') ~= 0 || exist('ratiodata','var') ~= 0
       if size(lengthdata,2) > 1 || size(widthdata,2) > 1 || size(ratiodata,2) > 1
           [centroid,outputs] = SOM_NOS(lengthdata,widthdata,ratiodata,angledata,nosdata,K);
           centromatrix(1,:,:) = centroid{1,1}; %centroid is a cell not a matix
       else
           continue
       end
   else 
       continue
   end
   
   
   %% result check
   % Identify the index of cluster
   for k = 1:size(outputs,1) %the number of class
       n =1;
       for j = 1:size(outputs,2)
           if outputs(k,j)~=0
               indexone(k,n) = j;
               n =n+1;
           end
       end
       hist(1,k) = n;
   end

    ccmatrix=zeros(K*numel(files),5); % 4 is number of feature
    n=1;
    for ii = 1:numel(files)
        for k = 1:K
            ccmatrix(n,:) = centromatrix(ii,k,:);
            n=n+1;
        end
    end
    
    % save the matrix into the excel file
    xlswrite(worksheetName,ccmatrix,'correlationMatrix');
    xlswrite(worksheetName,hist,'HistogramMeltpool');
    
    
    %% The statictical analysis ______%%%%%%%%%
    % K is the number of cluster
    for cluster = 1:K
        n=1;
        for i = 1:size(outputs,2)
            if outputs(cluster,i) ==1
                index(cluster,n) = i;
                n=n+1;
            end
        end
    end

    for cluster =1:K
        n=1;
        for in = 1:size(index(cluster,:),2)
            if index(cluster,in) ~= 0
                Length(cluster,n)= Data.DataBase(index(cluster,in)).Length;
                Width(cluster,n) = Data.DataBase(index(cluster,in)).Width;
                Ratio(cluster,n) = Data.DataBase(index(cluster,in)).Ratio;
                Angle(cluster,n) = Data.DataBase(index(cluster,in)).Angle;
                NOS(cluster,n) = Data.DataBase(index(cluster,in)).NOS;
                n=n+1;
            end
        end
    end
%% Cluster analysis(Find the abnormal and matrix)//Number of melt pool and correlatin matrix
% Using "Data" structure extract the data in order to analysis
for i = 1:K % cluster its statistical analysis
    Dimension(i,1)=mean(Length(i,:));
    Dimension(i,2)=mean(Width(i,:));
    Dimension(i,3)=mean(Ratio(i,:));
    Dimension(i,4)=mean(Angle(i,:));
    Dimension(i,5)=mean(NOS(i,:));
    Dimension(i,6)=std(Length(i,:));
    Dimension(i,7)=std(Width(i,:));
    Dimension(i,8)=std(Ratio(i,:));
    Dimension(i,9)=std(Angle(i,:));
    Dimension(i,10)=std(NOS(i,:));
end
    xlswrite(worksheetName,Dimension,'worksheet5');
    xlswrite(worksheetName,index,'index_of_cluster');
    %read the excel file
    [~,sheet_name]=xlsfinfo(worksheetName)
    for k=1:numel(sheet_name)
        data_sheet{k}=xlsread(worksheetName,sheet_name{k});
    end

    Meltpool_dimension = data_sheet{1,1}(1:size(data_sheet{1,1},1),1:size(data_sheet{1,1},2))'; %16x 1 Table
    %Number_of_Meltpool = data{1,1}(1:16,1:3); % 16 x1 Table

    corr_melt = corr(Meltpool_dimension);
    
    %figure,surf(corr_melt),colorbar   %/max(max(a(1:16,1:16)))
    
    for i = 1:size(corr_melt,1)
        med(i) = median(corr_melt(:,i));
    end
    %------------------CRITERIA-----------------------%
    % chose the criteria of the correlation matrix
    corr_index = [];
    n =1;
    for row = 1: size(corr_melt,1) % row of cluster
        for col = 1: size(corr_melt,2)  % col of cluster
            if corr_melt(row,col)< 1 
                corr_index(1,n) = row;
                corr_index(2,n) = col; % only extract one row, why i wrtie 2 row becuse it has to be seperated otherwise it can't be recognized!
                corr_index(3,n) = corr_melt(row,col);    
                n = n+1;
            end
        end
    end
    
    % criteria of std deviaion
    B = sort(Dimension(:,8)); % standard of ratio, 
    for num = 1: size(Dimension,1)
        if Dimension(num,6) == B(1) 
            std_num = num; %% This is good melt pool (find the cluster which has most of melt pool)
        end
    end
    % spatte angle chosen(criteria)
    Spatter_index = [];
    n =1;
    for L = 1: size(Data.DataBase,2)
        if (Data.DataBase(L).Angle <20) % the spatter angle 
            Spatter_index(n) = L; %% these image contain big spatter angel -> bad melt pool
            n = n+1;
        end
    end
    % Abnormal melt pool cluster
    abnormal_meltpool = [];
    n =1;
    if ~isempty(corr_index) && ~isempty(Spatter_index)
        for cluster_num = 1:size(corr_index,2)
            for melt_num = 1: size(index,2)
                for spatter_num = 1:size(Spatter_index,2)
                    if index(corr_index(1,cluster_num),melt_num) == Spatter_index(spatter_num) % if spatter numer in the melt pool number 
                        abnormal_meltpool(n,1) = index(cluster,melt_num); % melt pool number
                        abnormal_meltpool(n,2) = corr_index(3,cluster_num); % the correlation number of melt pool
                        n = n+1;
                    end
                end
            end
        end
    end
    % removeDuplicates
    if ~isempty(abnormal_meltpool)
        ab_index =1;
        for i = 1:size(abnormal_meltpool,1)
            if abnormal_meltpool(i,1) ~= abnormal_meltpool(ab_index,1)
                ab_index = ab_index +1;
                abnormal_meltpool(ab_index,1) = abnormal_meltpool(i,1);
            end
        end
    end

    %% Prediction Model combine the 3D result and validation (only one layer)
    for num_data = 1:size(Data.DataBase,2)
        if ndims(Data.DataBase(num_data).Image) == 3
            image = double(rgb2gray(Data.DataBase(num_data).Image));
        else
            image = double(Data.DataBase(num_data).Image);
        end
        [M_X,inde_x] = max(image);
        [m_x,ccentroid_x] = max(M_X);
    
        [M_Y,inde_y] = max(image');
        [m_y,ccentroid_y] = max(M_Y);
        ccentroid(1,num_data) = (ccentroid_x); %*13.67;% um/pixel
        ccentroid(2,num_data) = ccentroid_y;
        %%%%% no write in
        n =1;
        if isempty(abnormal_meltpool) == 0
            for NumberOfMeltpool = abnormal_meltpool(:,1)'
                if  num_data == NumberOfMeltpool
                    position(1,n,Layer) = ccentroid_x;
                    position(2,n,Layer) = ccentroid_y;
                    position(3,n,Layer) = abnormal_meltpool(find(num_data== NumberOfMeltpool),2);
                    n = n + 1;
                end
            end
        else
            position(1,1,Layer) = 0;
            position(2,1,Layer) = 0;
            position(3,1,Layer) = 0;
        end

    end
    % record the layer of position
    Layer = Layer + 1;
    
    clear('corr_melt');clear('outputs');clear('leng');clear('width');clear('ratio');clear('angle');clear('nos');
    clear('index');clear('lengthdata');clear('widthdata');clear('ratiodata');clear('angledata');clear('nosdata');
    clear('indexone');clear('hist');clear('Length');clear('Width');clear('Ratio');clear('NOS');
    clear('Angle');clear('Dimention');clear('data_sheet');clear('med');
    clear('corr_index');clear('Spatter_index');clear('abnormal_meltpool');clear('ccentroid');
end

%% extract the "position matrix to 3-D matrix"
dim_3D = [];
n = 0;
for i = 1: size(position,3)
    for ii = 1:size(position,2)
        dim_3D(n+ii,1) = position(1,ii,i);
        dim_3D(n+ii,2) = position(2,ii,i);
        dim_3D(n+ii,3) = i;
        dim_3D(n+ii,4) = position(3,ii,i);
    end
    n = n + ii;
end
%% Dara cleaning
dim_3D(dim_3D(:,2)== 0,:) = [];
dim_3D(dim_3D(:,1)== 0,:) = [];
% dim_3D(dim_3D(:,4) > 0.79,:) = [];
x = 0:size(Data.DataBase(1).Image,1); % row   
y = 0:size(Data.DataBase(1).Image,2); % col  
z = 1:numel(files); 

figure,scatter3(dim_3D(:,1),dim_3D(:,2),dim_3D(:,3),40,dim_3D(:,4),'filled');

%% density calculation
% Object volume (unit micro meter)
Long = 60;
Width = 2;
Height = 2;
Volume = Long * Width * Height;
% Pore volume summation
pore_volume = sum(dim_3D(:,4))*0.01;
 
Relative_density = abs(((Volume - pore_volume)/Volume))*100;



