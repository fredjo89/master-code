clc; clear all; close all;

nx = 28;
ny = 20;

G = cartGrid([nx, ny]);
G = computeGeometry(G);

% Create Coarse Grid
pv = partitionUI(G, [7, 4]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);


basis = zeros(nx*ny,1);

for i=1:28
    currentBasis = CG.cells.interaction{i}; 
    currentBoundary = findSupportBoundary(currentBasis,nx,ny);
    for j=1:nx*ny
       if currentBoundary(j)==1
           basis(j)=2;
       end
    end
    
 
    
end

for i=1:nx*ny
     if any(abs(CG.cells.centers-i)<1e-10)
        basis(i) = 4; 
     end
end




fileID = fopen('basisDistr.txt');
C = textscan(fileID,'%f');
fileData = C{1}; 
basisDistr = fileData(2:length(fileData))+1;

theDistr = zeros(nx*ny,1);
for i = 1:nx*ny
    basisNumber = pv(i);
    theDistr(i) = basisDistr(basisNumber);
end


colorMatrix = [ 1 1 1; .35 .35 .35; 1 .3 0; 0 .3 1; 0 0.75 0];
%FigHandle = figure('Position', [100, 400, 1200, 1200]);
figure();
hold on; 
plotCellData(G,basis,'EdgeColor', 'y')
plotGrid(G,'FaceColor', 'none')
%outlineCoarseGrid(G,pv,'k', 'lineWidth', 4)
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
%colormap(1,1) = 1

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters gridBoundary_partition


%{
%FigHandle = figure('Position', [1300, 400, 1200, 1200]);
figure();
plotCellData(G,theDistr,'EdgeColor', 'y')
plotGrid(G,'FaceColor', 'none')
%outlineCoarseGrid(G,pv,'k', 'lineWidth', 4)
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
caxis([.5 max(theDistr)+.5]);
%colormap(lines(max(theDistr)));
colormap(lines( length(unique(basisDistr)) +6   ));
%colormap(lines(28));


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters part_kway_28


disp(' ')
for i = 4 :-1: 1
    disp(basisDistr( (i-1)*7+1 : i*7)')
    
end
    

length(unique(basisDistr))
     %}
     
