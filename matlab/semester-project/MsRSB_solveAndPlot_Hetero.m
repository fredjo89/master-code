%% Hetero case
clc; clear all; close all; 

nBlocksX = 5;
nBlocksY = nBlocksX; 
nCellsX = 50; 
nCellsY = nCellsX; 
lengthX = nCellsX; 
lengthY = lengthX; 

iterations = 500; 

G = cartGrid([nCellsX, nCellsY], [lengthX, lengthY]);
G = computeGeometry(G);

p = gaussianField(G.cartDims, [.2 .4], [3 3 3] , 10.5);
K = p.^3*10^(-10)./(0.81*72*(1-p).^2);
rock.poro = p(:);
rock.perm = K(:);

gravity reset off; 

fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

pv = sum(poreVolume(G,rock));
src = addSource([], 1,-pv);
src = addSource(src, G.cells.num-nCellsX+1,pv); 
src = addSource(src, nCellsX,pv); 
src = addSource(src, G.cells.num, -pv); 

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

figure(); 
plotCellData(G,rock.poro);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 

figure(); 
plotCellData(G,rock.perm);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128)); 

error = state_fs.pressure - state_ms.pressure; 
error = abs(error); 

infNorm = max(error)/max(abs(state_fs.pressure))

twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))



