function [centroid,outputs] = SOM_NOS(lengthdata,widthdata,ratiodata,angledata,NOSdata,OutputDimension)
inputs = [lengthdata;widthdata;ratiodata;angledata;NOSdata];
%inputs=[irisInputs];

% Create a Self-Organizing Map
dimension1 = OutputDimension;
%dimension2 is output <worse good best>
dimension2 = 1;
distanceFcn ='linkdist';
coverSteps = 10;
initNeighbor = OutputDimension;
topologyFcn = 'hextop';
net = selforgmap([dimension1 dimension2],coverSteps,initNeighbor,topologyFcn,distanceFcn);
% Trianing parameter setup
net.trainParam.epochs = 1000;




% Train the Network
[net,tr] = train(net,inputs);

% Test the Network
outputs = net(inputs);
centroid = net.IW
% View the Network
%view(net)

%figure, plotsomtop(net)
%figure, plotsomhits(net,inputs)
%figure, plotsompos(net,inputs)
end
% Plots
% Uncomment these lines to enable various plots.
% figure, plotsomtop(net)
% figure, plotsomnc(net)
% figure, plotsomnd(net)
% figure, plotsomplanes(net)
% figure, plotsomhits(net,inputs)
% figure, plotsompos(net,inputs)