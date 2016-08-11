clc; clear all; close all;

nx = 27; 
ny = 27; 
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

CG = setupMexInteractionMapping(CG);
[offsets, support, celltypes] = getGridData(CG);

M = CG.cells.num; 
N = G.cells.num;

basis = zeros(1,N);
for i = offsets(5):offsets(6)
    if (celltypes(i)==1)
        basis(support(i)+1) = 1;
    else
        basis(support(i)+1) = 2;
    end
end

e = 0.5;
colorMatrix = [ 1 1 1; e e e ; 0.3  0.5  1    ];

hold on; 
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,basis','EdgeColor', 'y')
plotGrid(G,'FaceColor', 'none')
%outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1;
%set(colorbar, 'YTick',1:max(basis+1));

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters example-grid
%}


