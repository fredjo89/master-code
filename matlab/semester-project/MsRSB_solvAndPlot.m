%% Homo case
clc; clear all; close all; 

nBlocksX = 4;
nBlocksY = nBlocksX; 
nCellsX = 80; 
nCellsY = nCellsX; 
lengthX = nCellsX; 
lengthY = lengthX; 

iterations = 500; 


G = cartGrid([nCellsX, nCellsY], [lengthX, lengthY]);
G = computeGeometry(G);

p = ones(1, nCellsX*nCellsY);
K = ones(1, nCellsX*nCellsY);
rock.poro = p(:);
rock.perm = K(:);

gravity reset off; 

fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);


% Creating CG
pv = partitionUI(G,[nBlocksX,nBlocksY]); 
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);

poreVolume = 10*sum(poreVolume(G,rock));
src = addSource([], 1,-2*poreVolume);
src = addSource(src, G.cells.num-nCellsX+1,2*poreVolume); 
src = addSource(src, nCellsX,2*poreVolume); 
src = addSource(src, G.cells.num, -2*poreVolume); 

src = addSource(src, CG.cells.centers(nBlocksX*nBlocksY/2 +nBlocksX/2), -4*poreVolume); 
src = addSource(src, CG.cells.centers(nBlocksX*nBlocksY/2 -nBlocksX/2+1), 4*poreVolume); 

initState = initResSol(G,0.0, 1.0);

 

hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);



% Fine-scale solver
state_fs = incompTPFA(initState,G,hT, fluid, 'src', src); 


% MsRSB solver
basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',iterations);
state_ms = incompMultiscale(initState,CG,hT,fluid,basis,'src',src);


% Plotting results 
figure(); 
plotCellData(G,state_fs.pressure);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 

figure(); 
plotCellData(G,state_ms.pressure);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 

error = abs(state_fs.pressure - state_ms.pressure); 


infNorm = max(abs(error))/max(abs(state_fs.pressure))

twoNorm_fs = sqrt(sum(state_fs.pressure.^2))
twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig
%}





