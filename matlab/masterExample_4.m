clc; clear all; close all; 

% Number of cells
nx = 60;
ny = 30; 
nz = 10; 
% Number of coarse blocks
NX = 10;
NY = 5;
NZ = 2;

G = processGRDECL(simpleGrdecl([nx ny nz], 0.0));
G = computeGeometry(G);
%rock.perm = logNormLayers(G.cartDims, [100 400 50 350], 'indices', [1 2 5 7 11]);
rock.perm = logNormLayers(G.cartDims, [100 5000 ], 'indices', [1 6 11], 'std', 40.5, 'sigma', 15, 'a',0.1, 'sz', [9,3,3]);
%rock.perm = ones(G.cells.num,1);

%rock.perm = ones(G.cells.num,1);
gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

% Boundary conditions 
bc = pside([], G, 'left', 1);
bc = pside(bc, G, 'right', 0);

initState = initResSol(G,0.0, 1.0);
hT = computeTrans(G,rock); 

% Fine-scale solver
state_fs = incompTPFA(initState,G,hT, fluid,  'bc', bc); 

%rock.perm = ones(G.cells.num,1);
%hT = computeTrans(G,rock); 
%state_fs_2 = incompTPFA(initState,G,hT, fluid,  'bc', bc); 
%figure()
%plot((state_fs.pressure-state_fs_2.pressure)./state_fs.pressure);
%sum(abs((state_fs.pressure-state_fs_2.pressure)./state_fs.pressure))


% MsRSB solver
A = getIncomp1PhMatrix(G, hT);
%p = partitionUI(G,[NX NY NZ]); 
p = partitionUniformPadded(G, [NX NY NZ]);
%G_fault = makeInternalBoundary(G, find(G.faces.tag > 0))
%p = processPartition(G_fault, p);

CG = generateCoarseGrid(G,p); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegion(CG); 
CG = setupMexInteractionMapping(CG);

x = [ 0 5 10 20 50 75 100  ]; 
%x = 20;
for i = 1: length(x)

basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',x(i));
state_ms = incompMultiscale(initState,CG,hT,fluid,basis, 'bc',bc);

%state_ms.pressure = state_ms.pressure - min(state_ms.pressure); 
%state_ms.pressure= state_ms.pressure/max(state_ms.pressure);

% Calculating error norms
error = abs(state_fs.pressure - state_ms.pressure); 
infNorm(i) = max(error)/max(abs(state_fs.pressure));
twoNorm_error(i) = sqrt(sum(error.^2)/sum(state_fs.pressure.^2));
end

FigHandle = figure('Position', [1500, 200, 500, 500]);
plot(x,infNorm,'x--')
%axis([0,200,0,1]);
FigHandle = figure('Position', [2000, 200, 500, 500]);
plot(x,twoNorm_error,'x--')
%axis([0,200,0,1]);


%% Plots


FigHandle = figure('Position', [200, 1000, 600, 600]);
plotCellData (G , log10 (rock.perm ), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p, 'linewidth', 3)
colormap(jet(128)); 


FigHandle = figure('Position', [1400, 1000, 600, 600]);
plotCellData (G , full(basis.B(:,15)), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p, 'linewidth', 3)
colormap(jet(128)); 

FigHandle = figure('Position', [800, 1000, 600, 600]);
plotCellData (G , state_fs.pressure, 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p, 'linewidth', 3)
colormap(jet(128)); 


%{
FigHandle = figure('Position', [200, 200, 600, 600]);
plotCellData (G , abs(state_fs.pressure-state_ms.pressure), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p, 'linewidth', 3)
colormap(jet(128)); 
%}




FigHandle = figure('Position', [800, 200, 600, 600]);
plotCellData (G , state_ms.pressure, 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p)
colormap(jet(128)); 
outlineCoarseGrid(G,p, 'linewidth', 3)

%{

FigHandle = figure('Position', [200, 200, 1500, 1500]);
plotCellData (G , full(basis.B(:,1)), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p)
colormap(jet(128)); 
outlineCoarseGrid(G,p, 'linewidth', 3)



error = state_fs.pressure - state_ms.pressure; 
error = abs(error); 

infNorm = max(error)/max(abs(state_fs.pressure))

twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))
%}
