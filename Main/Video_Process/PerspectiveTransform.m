function tf=PerspectiveTransform(p)%% I is input img, p is corner coordinate
boardSize=[4,4];
squareSize = 160; % 轉換後方格紙邊長長度(pixel unit)
%imageSize = [size(I,1),size(I,2)];
%% define world points(the matrix that stores the target coordinates )
worldPoints = generateCheckerboardPoints(boardSize,squareSize);
%%test
% p=zeros(9,2);
% p(1,:)=[153 15];
% p(7,:)=[446 34];
% p(3,:)=[150 200];
% p(9,:)=[442 222];
% p(2,:)=(p(1,:)+p(3,:))/2;
% p(4,:)=(p(1,:)+p(7,:))/2;
% p(5,:)=(p(3,:)+p(7,:))/2;
% p(6,:)=(p(3,:)+p(9,:))/2;
% p(8,:)=(p(7,:)+p(9,:))/2;
%% define transformation matrix

tf = fitgeotrans(p,worldPoints,'projective');

% format shortg
% disp('tf.T');
% disp(tf.T)
% format
% fprintf('norm T: %g\n',norm(tf.T));
% fprintf('det T: %g\n',det(tf.T));

