clc; clear all; close all;

nx = 4*9; 
ny = 4*9; 
nCoarse = 4;
nBasis = 6; 
nBasis2 = 11; 


G = cartGrid([nx, ny]);
G = computeGeometry(G);


% Create Coarse Grid
pv = partitionUI(G, [nCoarse, nCoarse]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);


currentBasis = CG.cells.interaction{nBasis}; 
currentBasis2 = CG.cells.interaction{nBasis2};
for i=1:nx*ny
    if any(abs(currentBasis-i)<1e-10) && any(abs(currentBasis2-i)<1e-10)
        basis(i) = 1; 
    elseif any(abs(currentBasis-i)<1e-10)
        basis(i)=2; 
    elseif any(abs(currentBasis2-i)<1e-10)
        basis(i) = 3; 
    else
        basis(i) = 0; 
    end
    
     if any(abs(CG.cells.centers-i)<1e-10)
        basis(i) = 4; 
    end
    
end



%{
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
currentBasis = CG.cells.interaction{5}; 
currentBoundary = findSupportBoundary(currentBasis,nx,ny);
for i=1:nx*ny
    if currentBoundary(i)==1
        basis(i)=3;
    end
end


%}


my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;



%colorMatrix = [ 1 1 1; .55 .55 .55; my_red_2; my_blue_1; my_green_2];



colorMatrix = [ 1 1 1; .35 .35 .35; 1 .3 0; 0 .3 1; 0 0.75 0]
colorMatrix = [ 1 1 1; .35 .35 .35; my_red_3; my_blue_1; my_green_2]


hold on; 
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,basis','EdgeColor', 'k', 'EdgeAlpha', 0.01)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k', 'linewidth',3)
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1
%set(colorbar, 'YTick',1:max(basis+1));





