function [T,seg,result] = spatterSeg(I)
 seg=zeros(size(I,1),size(I,2));
for T=0.4:0.01:0.99
    BW1=imbinarize(I,T);%影像二值化
    [label,~]= bwlabel (BW1, 8);
    Label = regionprops (label,'Area','PixelList');  %%% mark calibration bar
    areaNum=size(Label,1); %紀錄分割的圖案個數
    minArea=min([Label.Area]);
    maxArea = max([Label.Area]);
    
    if size(Label,1)~=1
        [~,minidx] = min([Label.Area]);
        for a = 1:areaNum
             Areasmall = Label(minidx).PixelList;
            for i=1:size(Areasmall,1)
                seg(Areasmall(i,2),Areasmall(i,1))=255;
            end
        end
        break;
    end
end
seg=uint8(seg);

se=strel('disk',5);
seg=imdilate(seg,se);
result=I-seg;


