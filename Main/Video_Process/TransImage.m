% imagine Transform
function Imagine = TransImage(originalimagine,x1,y1,x3,y3,x7,y7,x9,y9)

p=zeros(9,2);
p(1,:)=[x1 y1];
p(3,:)=[x3 y3];
p(7,:)=[x7 y7];
p(9,:)=[x9 y9];
p(2,:)=(p(1,:)+p(3,:))/2;
p(4,:)=(p(1,:)+p(7,:))/2;
p(5,:)=(p(3,:)+p(7,:))/2;
p(6,:)=(p(3,:)+p(9,:))/2;
p(8,:)=(p(7,:)+p(9,:))/2;
tf=PerspectiveTransform(p);
   
[Imagine,~]=imwarp(originalimagine,tf);%原始影像座標轉換
end