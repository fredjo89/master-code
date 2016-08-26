clc; clear all;  

fileID = fopen('basisDistr.txt');
C = textscan(fileID,'%f');

fileData = C{1}; 

n = sqrt(fileData(1));

basisDistr = fileData(2:1+n^2);


G = cartGrid([6, 11], [2,1]);
G = computeGeometry(G);

figure(); 
hold on; 
plotCellData(G,basisDistr);
%plotCellData(G,redBlack,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
axis equal tight off; 
colormap(jet(1000)); 
outlineCoarseGrid(G,basisDistr,'k');
%colormap(colorMatrix)








