%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare the video 
%Tack read the video and ignore the black figure
%Create: 2021/02/04
%Editor: Castro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc,clear
close all

%% Perspective setting(cubic dimensionm and coordinate)

el = 5000;
x1 = 130;y1 = 201;
x3 = 137;y3 = 319;
x7 = 377;y7 = 199;
x9 = 386;y9 = 319;
pl = el/(x7-x1);

%% read the video
% folder = 'E:\0715EXP';
% path = dir(fullfile(folder,'1.avi'));

%---------------the main origin file work path-----------------------------------------
% fill the folder address in here which were save video file
%folder = 'F:\Matlab_analysis_data\0510_avi';
folder = 'E:\0715exp_avi';
path = dir(fullfile(folder,'*.avi'));
ncount = 1;
Numdata = 1;



%%%%%%%%%%%%%%Testing code region%%%%%%%%%
% video = 58
% name = [];
%     n = 1;
%     for num = 1: size(path(video).name,2)
%         if isnan(str2double(path(video).name(num))) == 0
%             name(n) = str2double(path(video).name(num));
%             n = n+1;
%         else
%             break;
%         end
%     end
%     name = sprintf('%1.f',name);
%     Inputname = strcat(name,'.mat');
%%%%%%%%%%%%%%%%%%%%%%%

for video = 1:numel(path)

filename = fullfile(path(video).folder,path(video).name)
% filename = fullfile(path(video).folder,'Output_18.avi')
%filename = 'F:\Code\ouput_video\Output_26.avi'
if ~isfile(filename)
    ncount = ncount+1;
%obj = VideoReader('H:\Castro\B_para\AVI_B\B1.avi');
else
obj = VideoReader(filename);
% if mod(ncount,5)==1
%     angle = 90;
% elseif mod(ncount,5)==2
%     angle = 135;
% elseif mod(ncount,5)==3
%     angle = 180;
% elseif mod(ncount,5)== 4
%     angle = 225;
% elseif mod(ncount,5) == 0
%     angle = 270;
% end

nFrames = obj.numberOfFrames;
vidHeight = obj.Height;
vidWidth = obj.Width;
n =1;
frame_grabed=nFrames;
I = ones(vidHeight,vidWidth,frame_grabed);
I = uint8(I);
%% Preallocate movie structure (initializet the object function)
mov(1:nFrames)= struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'));

for k = 1:nFrames
    mov(k).cdata = read(obj,k);
end
% image transform(perspective matrix)

%% detect the transforming Image size and collect the completly melt pool shape
zeroI = TransImage((mov(1).cdata),x1,y1,x3,y3,x7,y7,x9,y9);
TransHeight = size(zeroI,1);
TransWidth = size(zeroI,2);
n = 1;
for k = 1:nFrames
   if sum(sum(double(mov(k).cdata)))>10000
       n = n+1;
   end
end
% create the database and it saved the raw image, melt pool length, melt pool width, melt pool ratio,
% spatter angle and number of spatter(NOS)
if video == 1
    shape(1:n) = ...
        struct('Image',zeros(TransHeight,TransWidth,'uint8'),...
                'Length',0,...
                'Width',0,...
                'Ratio',0);
    DataBase(1:n) = ...
        struct('Image',zeros(TransHeight,TransWidth,'uint8'),...
                'Length',0,...
                'Width',0,...
                'Ratio',0,...
                'Angle',0,...
                'NOS',0); % NOS is Number Of Spatter
else
    video;
end

n = 1;
%% delete the black image based on the sum of intensity
for k = 1:nFrames
   if sum(sum(double(mov(k).cdata)))>10000
       shape(n).Image = mov(k).cdata;
       n = n+1;
   end
       
end

% analysis each image from "shape"

for each = 1: size(shape,2)

    if ndims(shape(each).Image)==3 
        Iold = rgb2gray(shape(each).Image);
    else
        Iold = shape(each).Image;
    end
    
    if sum(sum(double(Iold))) <= 21000 %% in ordr to negelct the weird melt pool dimension 
        continue;
    end
    
    TransfImage = TransImage(Iold,x1,y1,x3,y3,x7,y7,x9,y9);

    Tr = TransfImage;
    
    Binary_Iold = Iold;
    T = double(Tr);
    
    % Find the maximum algorithm
    [M,index] = max(T');
    [m,in] = max(M);
    TS = T(in,:); %Target Signal
    
    % filter using the 1-D Fourier Lowpass Filter
    L = fourierLowPass(TS,50,100);
    DL1 = diff(L);
    DL2 = diff(DL1);
    
    %Data Length
    xL1 = [1:1:size(L,2)]*(pl);
    xL2 = [1:1:size(DL1,2)]*(pl);
    xL3 = [1:1:size(DL2,2)]*(pl);
    
    %% find the maximum location of 1st derivative and watch the video
    [~,location] = max(DL1);
    threshold = L(1,location);
    
    for i = 1:size(Iold,1)
        for j = 1:size(Iold,2)
            if Iold(i,j) < threshold
                Binary_Iold(i,j) = 0;
            end
        end
    end
    
    for i = 1:size(Binary_Iold,1)
        for j = 1:size(Binary_Iold,2)
            if Binary_Iold(i,j)>0
                Binary_Iold(i,j)=255;
            end
        end
    end
    
    Binary_Iold = double(Binary_Iold);
     
    [label, number] = bwlabel(Binary_Iold, 8);
    Label = regionprops(label,'Area','Centroid','BoundingBox','PixelList');
    [~,idx] = max([Label.Area]);
    if isempty(idx) == 1
        each;
    else
        [G] = Label(idx).BoundingBox;
        length = G(3).*pl;
        width = G(4).*pl;
        
        final = Binary_Iold;
        seg = zeros(size(Binary_Iold,1),size(Binary_Iold,2));
        if size(Label,1)~=1
            for a = 1:number
                if Label(a).Area < Label(idx).Area
                    onesmall = Label(a).PixelList;
                    for i = 1:size(onesmall,1)
                        seg(onesmall(i,2),onesmall(i,1)) =255;
                    end
                    seg = uint8(seg);
                    se = strel('disk',5);
                    seg = imdilate(seg,se);
                    final = final - double(seg);
                end
            end
        end
    end

  if size(Label,1)> 2
      sort_Area = sort([Label.Area],'descend');
      for nn = 1: size(Label,1)
         if  Label(nn).Area == sort_Area(1) % the bigger one (melt pool)
             No_meltpool = nn;
         end
%          if Label(nn).Area == sort_Area(2) % the spatter
%              No_spatter = nn;
%          end
      end
      % Calculate spatter angle base upon the arctan function
      q =1;
      for Number = 1 : size(sort_Area,2) % Number is spatter number
          if Number ~= No_meltpool
            theta_matrix(q) = abs(atan((Label(Number).Centroid(2)-Label(No_meltpool).Centroid(2))/(Label(Number).Centroid(1)-Label(No_meltpool).Centroid(1)))*(180/pi));
            q = q + 1; % Number of spatter
          else
              continue;
          end
      end
      
      theta = mean(theta_matrix);
      theta_std = std(theta_matrix);
      
%       if isnan(theta) == 1
%           theta = 90;
%           theta_std = 0;
%       else
%         theta_std = std(theta_matrix);
%       end
      
  else
      theta = 0;
      q = 0;
  end
    



    DataBase( Numdata).Image = shape(each).Image;
    DataBase( Numdata).Length  = length;
    DataBase( Numdata).Width = width;
    DataBase( Numdata).Shape = final;
    DataBase( Numdata).Ratio = width/length;
    DataBase( Numdata).Angle = theta;
    DataBase( Numdata).NOS = q;

    clear Iold
    clear Binary_Iold
    clear T
    
    %%spatter segamentation?????
    %%U-net model or interpolation???
    Numdata = Numdata+1;
end
    name = [];
    n = 1;
    for num = 1: size(path(video).name,2)
        if isnan(str2double(path(video).name(num))) == 0
            name(n) = str2double(path(video).name(num));
            n = n+1;
        else
            break;
        end
    end
    % the amount of data name
    name = sprintf('%1.f',name);
    InputName = strcat(name,'.mat');
    
    
    %% testing file name
    %InputName = sprintf('test%d.mat',video);
    %InputName = name
%     InputName ='testing_4.mat';
    
    
    
    save(InputName,'DataBase')
    close all
    clear shape
    ncount=ncount+1;
end
end
% Avg_Length = mean([DataBase(1:end).Length])
% Avg_Width = mean([DataBase(1:end).Width])
% STD_Length = std([DataBase(1:end).Length])
% STD_Width = std([DataBase(1:end).Width])
