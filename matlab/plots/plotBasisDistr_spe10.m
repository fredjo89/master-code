clc; clear all; close all;


[G, rock, p, testcase] = makeProblem(); 
G = computeGeometry(G);

plotCellData(G,p);
outlineCoarseGrid(G,p);
view(3)

clear all; 

nx = 6;
ny = 11;
nz = 17; 

G = cartGrid([nx, ny, nz], [1200,2200,800]);
G = computeGeometry(G);




fileID = fopen('basisDistr.txt');
C = textscan(fileID,'%f');
fileData = C{1}; 
basisDistr = fileData(2:length(fileData))+1;

%{
theDistr = zeros(G.cells.num,1);
for i = 1:G.cells.num
    basisNumber = p(i);
    theDistr(i) = basisDistr(basisNumber);
end
%}




%FigHandle = figure('Position', [200, 200, 13*29, 11.5*29]);
figure();
%plotGrid(G,'FaceColor', 'none')
plotCellData(G,basisDistr)
%outlineCoarseGrid(G,p)
view(3)
caxis([.5 max(basisDistr)+.5]);
colormap(lines(max(basisDistr)));
axis equal tight off; 
%colormap(lines( length(unique(basisDistr)) +3   ));
%outlineCoarseGrid(G,pv,'k', 'lineWidth', 4)
%outlineCoarseGrid(G,p)
%view(3)
%caxis([.5 max(theDistr)+.5]);
%colormap(lines(max(theDistr)));
%colormap(lines( length(unique(basisDistr)) +6   ));
%colormap(lines(28));

%{
FigHandle = figure('Position', [200, 200, 13*150, 11.5*150]);
slice(reshape(theDistr, 60,220,10),[1 220],60,[1,85])
shading flat; axis equal off; set(gca,'zdir','reverse'), box on;
colorbar('horiz')
%}

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters distr_spe_20
%}


length(unique(basisDistr))
     
     
