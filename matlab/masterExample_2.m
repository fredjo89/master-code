clc; clear all; close all; 

% Number of coarse blocks (Upscaling factor is 125 atm)
NX = 10;
NY = 6;
NZ = 2;

G = processGRDECL(simpleGrdecl([50 30 10], 0.12));
G = computeGeometry(G);
%rock.perm = logNormLayers(G.cartDims, [100 400 50 350], 'indices', [1 2 5 7 11]);
rock.perm = logNormLayers(G.cartDims, [1 1 1 1], 'indices', [1 2 5 7 11], 'std', 100);
rock.perm = ones(G.cells.num,1);
gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

% Boundary conditions 
bc = pside([], G, 'left', 1);
bc = pside(bc, G, 'right', 0);

% Source terms (not used atm)
pv = G.cells.num;
src = addSource([], G.cells.num,-pv);
src = addSource(src, 1, pv); 
for i = 1:29
    src = addSource(src, G.cells.num-50*i, -pv); 
    src = addSource(src,50*i+1, pv);
end

initState = initResSol(G,0.0, 1.0);
hT = computeTrans(G,rock); 

% Fine-scale solver
state_fs = incompTPFA(initState,G,hT, fluid,  'bc', bc); 

% MsRSB solver
A = getIncomp1PhMatrix(G, hT);
p = partitionUI(G,[NX NY NZ]); 
%p = partitionUniformPadded(G, [NX NY NZ]);

G_fault = makeInternalBoundary(G, find(G.faces.tag > 0))
p = processPartition(G_fault, p);

CG = generateCoarseGrid(G,p); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG); 

basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',10);
state_ms = incompMultiscale(initState,CG,hT,fluid,basis, 'bc',bc);

% Calculating error norms
error = abs(state_fs.pressure - state_ms.pressure); 
infNorm = max(error)/max(abs(state_fs.pressure))
twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))



%% Plots
FigHandle = figure('Position', [200, 200, 1500, 1500]);
plotCellData (G , log10 (rock.perm ), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
plotGrid(G,src.cell, 'FaceColor', 'w'); 
outlineCoarseGrid(G,p, 'linewidth', 3)
colormap(jet(128)); 


FigHandle = figure('Position', [200, 200, 1500, 1500]);
plotCellData (G , state_fs.pressure, 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
plotGrid(G,src.cell, 'FaceColor', 'w'); 
colormap(jet(128)); 

%{
FigHandle = figure('Position', [200, 200, 1500, 1500]);
plotCellData (G , state_ms.pressure, 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
plotGrid(G,src.cell, 'FaceColor', 'w'); 
outlineCoarseGrid(G,p)
colormap(jet(128)); 
outlineCoarseGrid(G,p, 'linewidth', 3)


FigHandle = figure('Position', [200, 200, 1500, 1500]);
plotCellData (G , full(basis.B(:,1)), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
plotGrid(G,src.cell, 'FaceColor', 'w'); 
outlineCoarseGrid(G,p)
colormap(jet(128)); 
outlineCoarseGrid(G,p, 'linewidth', 3)



error = state_fs.pressure - state_ms.pressure; 
error = abs(error); 

infNorm = max(error)/max(abs(state_fs.pressure))

twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))
%}
