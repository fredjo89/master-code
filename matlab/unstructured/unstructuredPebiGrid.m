mrstModule add voronoi2d

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

%% Alternative 2: Use cart partitioning, and subsample
% We create a cartesian grid of same size, and then sample the partition
% vector onto the new grid

%G_cart = cartGrid(k_dim);
%p = partitionUI(G_cart, [3, 3, 2]);
%p = sampleFromBox(G, reshape(p, k_dim));

figure; plotCellData(G, p)
view(-40, 50)

