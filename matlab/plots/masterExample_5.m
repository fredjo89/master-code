clc; clear all; close all; 

% Number of cells
nx = 60;
ny = 30; 
nz = 10; 
% Number of coarse blocks
NX = 5;
NY = 3;
NZ = 2;

G = processGRDECL(simpleGrdecl([nx ny nz], 0.0));
G = computeGeometry(G);
%rock.perm = logNormLayers(G.cartDims, [100 400 50 350], 'indices', [1 2 5 7 11]);

%rock.perm = ones(G.cells.num,1);

makeRock = 1; 

if makeRock==0
    RR = 6; 
    rock.perm = logNormLayers(G.cartDims, [1000 ], 'indices', [1 11], 'std', 4.5, 'sigma',17, 'a',0.6, 'sz', [15+RR,7+RR,5+RR]);
    fid=fopen('MyRock.txt','w');
    fprintf(fid, '%f \n', rock.perm);
    fclose(fid); 
else
    fileID = fopen('MyRock.txt');
    C = textscan(fileID,'%f');
    rock.perm = C{1}; 
end
%rock.perm = convertFrom(rock.perm,milli*darcy);





max(rock.perm)/min(rock.perm)

max(rock.perm)
min(rock.perm)

gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

% Boundary conditions 
bc = pside([], G, 'left', 1);
bc = pside(bc, G, 'right', 0);

initState = initResSol(G,0.0, 1.0);
hT = computeTrans(G,rock); 

% Fine-scale solver
state_fs = incompTPFA(initState,G,hT, fluid,  'bc', bc); 

% MsRSB solver
A = getIncomp1PhMatrix(G, hT);
%p = partitionUI(G,[NX NY NZ]); 
p = partitionUniformPadded(G, [NX NY NZ]);
G_fault = makeInternalBoundary(G, find(G.faces.tag > 0))
p = processPartition(G_fault, p);

CG = generateCoarseGrid(G,p); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegion(CG); 
CG = setupMexInteractionMapping(CG);

x = [ 0 2 4 7 10 15 20 25 30 40 50 60 70 80 90 100]; 
%x = 0;
for i = 1: length(x)
basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',x(i));
state_ms = incompMultiscale(initState,CG,hT,fluid,basis, 'bc',bc);

% Calculating error norms
error = abs(state_fs.pressure - state_ms.pressure); 
infNorm(i) = max(error)/max(abs(state_fs.pressure));
twoNorm_error(i) = sqrt(sum(error.^2)/sum(state_fs.pressure.^2));
end

%{
rock.perm = ones(G.cells.num,1);
hT = computeTrans(G,rock); 
state_fs_2 = incompTPFA(initState,G,hT, fluid,  'bc', bc); 
error = abs(state_fs.pressure - state_fs_2.pressure); 
max(error)/max(abs(state_fs.pressure))
sqrt(sum(error.^2)/sum(state_fs.pressure.^2))
%}





%% Plots

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;

FigHandle = figure('Position', [1200, 200, 13*29, 11.5*29]);
plot(x,infNorm,'--o','Color',my_green_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_green_1, 'MarkerFaceColor', my_green_1);
xlabel('Iterations');
ylabel('Error');
%LEG1 = legend('inf Norm');
axis([0,max(x),0,0.35]);
set(gca,'fontsize',15)

FigHandle = figure('Position', [1800, 200, 13*29, 11.5*29]);
plot(x,twoNorm_error,'--o','Color',my_blue_1, 'LineWidth', 1, ...
   'MarkerSize', 7, 'MarkerEdgeColor', my_blue_1, 'MarkerFaceColor', my_blue_1);
xlabel('Iterations');
ylabel('Error');
axis([0,max(x),0,0.20]);
%LEG1 = legend('Two Norm');
set(gca,'fontsize',15)








% Plotting permeability
FigHandle = figure('Position', [200, 1000, 1000, 600]);
%FigHandle = figure('Position', [1800, 200, 13*29, 11.5*29]);
plotCellData (G , log10 (rock.perm ), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz');

ax = gca; 
axpos = ax.Position;
cpos = h.Position;

cpos(1) = 0.27;
cpos(3) = 0.65*cpos(3);
cpos(4) = 1.5*cpos(4);

h.Position = cpos;
ax.Position = axpos;

ticks = [  100 1000 10000]
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p, 'linewidth', 4)
colormap(jet(128)); 
set(gca, 'FontSize',55);



%{
% Plotting a basis function
FigHandle = figure('Position', [1400, 1000, 600, 600]);
plotCellData (G , full(basis.B(:,9)), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p, 'linewidth', 3)
colormap(jet(128)); 
%}

% Plotting fine scale pressure solution
FigHandle = figure('Position', [1000, 1000, 1000, 600]);
plotCellData(G , state_fs.pressure, 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz');
ax = gca; 
axpos = ax.Position;
cpos = h.Position;
cpos(1) = 0.27;
cpos(3) = 0.65*cpos(3);
cpos(4) = 1.5*cpos(4);
h.Position = cpos;
ax.Position = axpos;
ticks =[ 0.1 0.5 0.9 ];
set (h , 'XTick',  ticks , 'XTickLabel',ticks );
colormap(jet(128)); 
set(gca, 'FontSize',55);













%{
% Plotting difference between multiscale and finescale pressure
FigHandle = figure('Position', [200, 200, 600, 600]);
plotCellData (G , abs(state_fs.pressure-state_ms.pressure), 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p, 'linewidth', 3)
colormap(jet(128)); 



% Plotting multiscale pressure solution
FigHandle = figure('Position', [800, 200, 1000, 600]);
plotCellData (G , state_ms.pressure, 'EdgeColor','k'); view (45,30);
axis tight off , set ( gca , 'DataAspect',[0.5 1 1])
h= colorbar ('horiz'); ticks =25*2.^[0:5];
set (h , 'XTick', log10 ( ticks ), 'XTickLabel',ticks );
outlineCoarseGrid(G,p)
colormap(jet(128)); 
outlineCoarseGrid(G,p, 'linewidth', 3)
%}


%}

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters TRY_Johansen
%}


%
%print -dpng -r1000 TEMPFIG2

%print -dpng unstructuredGrid



