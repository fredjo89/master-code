%% Hetero case
clc; clear all; close all; 

nBlocksX = 10;
nBlocksY = nBlocksX; 
nCellsX = 100; 
nCellsY = nCellsX; 
lengthX = nCellsX; 
lengthY = lengthX; 

%lognormlayers

iterations = 500; 

G = cartGrid([nCellsX, nCellsY], [lengthX, lengthY]);
G = computeGeometry(G);


rock.perm = logNormLayers([nCellsX, nCellsY])
rock.poro = rock.perm / max(rock.perm)


gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

pv = sum(poreVolume(G,rock));
pv = 1; 

src = addSource([], 1,-pv);
src = addSource(src, G.cells.num, -pv); 
for (i = 1:nCellsY-1)
    src = addSource(src, i*nCellsX+1,-pv); 
    src = addSource(src, i*nCellsX,pv); 
end

 


initState = initResSol(G,0.0, 1.0);

hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);

% Fine-scale solver
state_fs = incompTPFA(initState,G,hT, fluid, 'src', src); 

% Creating CG
pv = partitionUI(G,[nBlocksX,nBlocksY]); 
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG); 

% MsRSB solver
basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',iterations);
state_ms = incompMultiscale(initState,CG,hT,fluid,basis,'src',src);


% Plotting results 
FigHandle = figure('Position', [200, 1200, 500, 500]);
plotCellData(G,state_fs.pressure);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 

FigHandle = figure('Position', [1000, 1200, 500, 500]);
plotCellData(G,state_ms.pressure);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 

FigHandle = figure('Position', [200, 200, 500, 500]);
plotCellData(G,rock.poro);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 

FigHandle = figure('Position', [1000, 200, 500, 500]);
plotCellData(G,rock.perm);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 


error = state_fs.pressure - state_ms.pressure; 
error = abs(error); 

infNorm = max(error)/max(abs(state_fs.pressure))

twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))




















