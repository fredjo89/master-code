close all; 

mrstModule add voronoi2d;

% Dimensions of grid - make box of 1000 by 1000 meter
pdim = 1000;
% zdim is the vertical size of the grid
zdim = 300;
% Number of cells in x/y direction
n_xy = 4*3;
% Number of layers
n_z = n_xy/4;
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




newP = zeros(G.cells.num, 1); 



for (i = 1:G.cells.num)
   newP(i) = randi([1 G.cells.num],1,1);
end


% Plot final result
figure;
plotCellData(G, newP)
plotGrid(G,'FaceColor','none','linewidth',2);
axis equal tight off; 
%colormap(jet(16));
view(-40, 30)














