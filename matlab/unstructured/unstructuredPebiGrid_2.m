close all; 

mrstModule add voronoi2d;

% Dimensions of grid - make box of 1000 by 1000 meter
pdim = 1000;
% zdim is the vertical size of the grid
zdim = 300;
% Number of cells in x/y direction
n_xy = 4*5;
% Number of layers
n_z = 4;
% Compute desired cell size
resGridSize = pdim./n_xy;

G_2d = pebiGrid(resGridSize, [pdim, pdim]);
% Extrude 2D unstructured grid to 3D semi-unstructured
G = makeLayeredGrid(G_2d, n_z);
% Scale nodes to get correct dimensions in z-direction
G.nodes.coords(:, 3) = zdim*G.nodes.coords(:, 3)./n_z;
G = computeGeometry(G);
%% Plot grid
figure;
plotGrid(G)
view(-40, 50)

%% Create perm and sample to grid
k_dim = [n_xy, n_xy, n_z];

% Define lognormal permeability, with means of 50, 20 and 100
rng(0);
K = logNormLayers(k_dim, [50, 20, 100]);
k = sampleFromBox(G, reshape(K, k_dim));
% Multiply by millidarcy and assign to rock
rock = makeRock(G, k*milli*darcy, 0.3);

% Plot final result
figure;
plotCellData(G, log10(rock.perm))
view(-40, 50)
%% First option: Use METIS
n_c = 5;
T = computeTrans(G, rock);
mrstModule add coarsegrid incomp

% First option: Use Metis if installed. Need to define METISPATH
% if METIS is found under /usr/bin
global METISPATH
METISPATH = '/usr/local/bin/';
%p = partitionMETIS(G, T, n_c);
p = ones(G.cells.num,1);


% My code

A = getIncomp1PhMatrix(G, T);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG);
CG = setupMexInteractionMapping(CG);

fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', 'unstructured']);
%fn = fullfile('/global/work/fredjoha/mrst_data', ['basis_', lower(testcase)]);


temp = zeros(CG.cells.num, 1); 
newP = zeros(G.cells.num, 1); 
for i = 1:CG.cells.num
    temp(i) = randi([1 CG.cells.num],1,1);
end


for (i = 1:G.cells.num)
   newP(i) = temp(p(i));
end

%{
FigHandle = figure('Position', [1200, 200, 800, 800]);
plotCellData(G,p,'EdgeColor', 'k','EdgeAlpha', 0.1)
view(-40, 50)
outlineCoarseGrid(G,p,'linewidth',2)
axis tight off; 
%}

%{
FigHandle = figure('Position', [1200, 200, 800, 800/(1.5)]);
plotCellData(G,newP,'EdgeColor', 'k','EdgeAlpha', 0.1)
view(-40, 55)
outlineCoarseGrid(G,p,'linewidth',2)
axis tight off; 
%}


%{
figure();
%view(-15,40);
%plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotCellData(G,boundary, find(boundary>0.1), ...
             'EdgeColor','k','EdgeAlpha',0.1);
colormap(colorMatrix)
colormap(1,1) = 1
%}

%print -dpng 3D_example_125









%cppMultiscaleBasis(CG, A, 'doSolve', false, 'writePath', fn);
newWrite2(CG, A, 'writePath', fn);


disp('Stored files:')
ls(fullfile(fn, 'input'))








