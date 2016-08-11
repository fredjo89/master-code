clc; clear all; close all;

nx = 30; 
ny = 30; 
nCoarse = 3;
nBasis = 5; 
nBasis2 = 1; 


G = cartGrid([nx, ny]);
G = computeGeometry(G);


% Create Coarse Grid
pv = partitionUI(G, [nCoarse, nCoarse]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
%CG = storeInteractionRegionCart(CG);



currentBasis = CG.cells.interaction{nBasis}; 
currentBasis2 = CG.cells.interaction{nBasis2};
for i=1:nx*ny
    if any(abs(currentBasis-i)<1e-10)
        basis(i) = 1; 
    elseif any(abs(currentBasis-i)<1e-10)
        basis(i)=0; 
    elseif any(abs(currentBasis2-i)<1e-10)
        basis(i) = 0; 
    else
        basis(i) = 0; 
    end
    
    
    
    
end




%close all; 
globalBoundary = zeros(1,nx*ny);
%find global boundary cells
for i=1:nCoarse*nCoarse
    currentBasis = CG.cells.interaction{i}; 
    currentBoundary = findSupportBoundary(currentBasis,nx,ny);
    for i=1:nx*ny
       if currentBoundary(i)==1
           basis(i)=2;
       end
    end
end
currentBasis = CG.cells.interaction{nBasis}; 
currentBoundary = findSupportBoundary(currentBasis,nx,ny);
for i=1:nx*ny
    if currentBoundary(i)==1
        basis(i)=3;
    end
    
     if any(abs(CG.cells.centers-i)<1e-10)
        basis(i) = 4; 
    end
end

colorMatrix = [ 1 1 1; .55 .55 .55; 1 .3 0; 0 .3 1; 0 0.75 0];




hold on; 
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,basis','EdgeColor', 'y')
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1
%set(colorbar, 'YTick',1:max(basis+1));

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters tempFig
%}


