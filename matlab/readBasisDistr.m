clc; clear all;  

fileID = fopen('basisDistr.txt');
C = textscan(fileID,'%f');

fileData = C{1}; 

n = sqrt(fileData(1));

basisDistr = fileData(2:1+n^2);

G = cartGrid([n, n]);
G = computeGeometry(G);


myBlue = [ 0.1 0.1 0.8];
myWhite = [0.9 0.9 0.9];
myGray = [0.55 0.55 0.55];
%colorMatrix = [1 0 0; 0 0 0  ];

figure(); 
hold on; 
plotCellData(G,basisDistr);
%plotCellData(G,redBlack,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
axis equal tight off; 
colormap(jet(1000)); 
outlineCoarseGrid(G,basisDistr,'k');
%colormap(colorMatrix)









