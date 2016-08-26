close all; 

mrstModule add voronoi2d;

% Dimensions of grid - make box of 1000 by 1000 meter
pdim = 1000;
% zdim is the vertical size of the grid
zdim = 300;
% Number of cells in x/y direction
n_xy = 20;
% Number of layers
n_z = 5;
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
n_c = 25;
T = computeTrans(G, rock);
mrstModule add coarsegrid incomp

% First option: Use Metis if installed. Need to define METISPATH
% if METIS is found under /usr/bin
global METISPATH
METISPATH = '/usr/local/bin/';
p = partitionMETIS(G, T, n_c);

figure; plotCellData(G, p)
view(-40, 50)


% My code

A = getIncomp1PhMatrix(G, T);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegion(CG);
CG = setupMexInteractionMapping(CG);

fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', 'unstructured']);
%fn = fullfile('/global/work/fredjoha/mrst_data', ['basis_', lower(testcase)]);



%cppMultiscaleBasis(CG, A, 'doSolve', false, 'writePath', fn);
newWrite2(CG, A, 'writePath', fn);


disp('Stored files:')
ls(fullfile(fn, 'input'))









