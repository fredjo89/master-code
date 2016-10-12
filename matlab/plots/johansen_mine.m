clc; clear all; close all; 


dpath = getDatasetPath('johansen');
sector = fullfile(dpath, 'NPD5');
filename = [sector, '.grdecl'];

grdecl        = readGRDECL(filename);  clear filename
actnum        = grdecl.ACTNUM;
grdecl.ACTNUM = ones(prod(grdecl.cartDims),1);
G             = processGRDECL(grdecl, 'checkgrid', false);

%%
% Plot the results
clf, subplot('position',[0.025 0.025 0.95 0.95]);
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotFaces(G,find(G.faces.tag>0),'FaceColor','r');
axis tight off; view(-145,60);

%%
% Next we mark the active part of the model
plotGrid(G,find(actnum(G.cells.indexMap)), ...
         'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeAlpha', 0.1);
view(20,75);

%% Height map
% It is only meaningful to show a height map of the active cells.
% Therefore, to inspect only the active model, we reset the ACTNUM field to
% its original values and recreate the grid. Now, inactive cells will be
% ignored and we therefore get a different unstructured grid.
grdecl.ACTNUM = actnum; clear actnum;
G = processGRDECL(grdecl); clear grdecl;
G = computeGeometry(G);

% Plotting a height map of the field using the z-component of the centroids
% of the cells
clf,
plotCellData(G,G.cells.centroids(:,3),'EdgeColor','k','EdgeAlpha',0.1);
colorbar, view(3), axis tight off, view(-20,40), zoom(1.2)

%% Porosity
% The porosity data are given with one value for each cell in the model. We
% read all values and then pick only the values corresponding to active
% cells in the model.
clf
p = reshape(load([sector, '_Porosity.txt'])', prod(G.cartDims), []);
poro = p(G.cells.indexMap); clear p
colorbar; caxis([0.1 0.3]), view(-45,15), axis tight off, zoom(1.2)



%% Permeability
% The permeability is given as a scalar field (Kx) similarly as the
% porosity. The tensor is given as K = diag(Kx, Kx, 0.1Kx) and we therefore
% only plot the x-component, Kx, using a logarithmic color scale.
clf
K = reshape(load([sector, '_Permeability.txt']')', prod(G.cartDims), []);
perm = bsxfun(@times, [1 1 0.1], K(G.cells.indexMap)).*milli*darcy; clear K;
rock = makeRock(G, perm, poro);




%%
% To show more of the permeability structure, we strip away the shale
% layers, starting with the layers with lowest permeability on top.
FigHandle = figure('Position', [700, 200, 1300, 1000]);
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
hp = plotCellData(G,log10(rock.perm(:,1)), ...
                  find(rock.perm(:,1)>0.01*milli*darcy), ...
                  'EdgeColor','k', 'EdgeAlpha', 1);
view(-20,35), axis tight off, zoom(1.2)      

% Manipulate the colorbar to get the ticks we want
h = colorbar;
cs = [0.01 0.1 1 10 100 1000];
caxis(log10([min(cs) max(cs)]*milli*darcy));
set(h, 'XTick', 0.8, 'XTickLabel','mD', ...
   'YTick', log10(cs*milli*darcy), 'YTickLabel', num2str(cs'));
set(gca,'fontsize',30)

colormap(jet(64)); 

ax = gca;
axpos = ax.Position;
cpos = h.Position;

cpos(3) = 2*cpos(3);
h.Position = cpos;
ax.Position = axpos;

x=get(h,'Position');
x(2)=0.22;
x(4)=0.65;
set(h,'Position',x)

set(get(h,'title'),'string','mD');



print -dpng TRY_Johansen
