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

dim_3D(dim_3D(:,2)== 0,:) = [];
dim_3D(dim_3D(:,1)== 0,:) = [];

x = 0:size(Data.DataBase(1).Image,1); % row   
y = 0:size(Data.DataBase(1).Image,2); % col  
z = 1:numel(files); 

figure,scatter3(dim_3D(:,1),dim_3D(:,2),dim_3D(:,3),40,dim_3D(:,4),'filled');
