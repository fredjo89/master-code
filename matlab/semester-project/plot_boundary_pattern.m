clc; clear all; close all;

nx = 16*4; 
ny = 9*4; 
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

% Create both support regions
currentBasis = CG.cells.interaction{nBasis}; 
currentBasis2 = CG.cells.interaction{nBasis2};
currentBasis3 = CG.cells.interaction{3};
currentBasis4 = CG.cells.interaction{7};
currentBasis5 = CG.cells.interaction{9};
for i=1:nx*ny
    if any(abs(currentBasis-i)<1e-10)
        basis(i) = 1; 
    else
        basis(i) = 0; 
    end
    if any(abs(currentBasis2-i)<1e-10)
        basis2(i) = 1; 
    else
        basis2(i) = 0; 
    end
    if any(abs(currentBasis3-i)<1e-10)
        basis3(i) = 1; 
    else
        basis3(i) = 0; 
    end
    if any(abs(currentBasis4-i)<1e-10)
        basis4(i) = 1; 
    else
        basis4(i) = 0; 
    end
    if any(abs(currentBasis5-i)<1e-10)
        basis5(i) = 1; 
    else
        basis5(i) = 0; 
    end
end


for i=1:nCoarse*nCoarse
    currentBoundary = findSupportBoundary(CG.cells.interaction{i},nx,ny);
    for i=1:nx*ny
       if currentBoundary(i)==1 && basis(i)==1
          basis(i)=2;
       end
    end
end

for i=1:nx*ny
    if basis(i)==1
        basis(i)=0;
    elseif basis(i) ==2
        basis(i) = 1; 
    end
end


for i=1:nx*ny
    if basis2(i)==1 && basis(i)~=2
        basis(i)=2;
    end
    if basis3(i)==1 && basis(i)~=2
        basis(i)=5;
    end
    if basis4(i)==1 && basis(i)~=2
        basis(i)=4;
    end
    if basis5(i)==1 && basis(i)~=2
        basis(i)=3;
    end

end




colorMatrix = [ 1 1 1 ; 1 0 0; ...
            0 0.4 0; 0. 1 0; ...
            0 0 1; 0 0 0.5];



hold on; 
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,basis','EdgeColor', 'y')
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1
%set(colorbar, 'YTick',1:max(basis+1));





