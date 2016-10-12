clc; clear all; close all;



%boxSize = 5; 
%n_box = 4; 

boxSize = 10; 
n_box = 2; 






G = cartGrid([boxSize*n_box, boxSize*n_box, boxSize*n_box]);
G = computeGeometry(G);

% Create Coarse Grid
pv = partitionUI(G, [n_box, n_box n_box]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
%CG = storeInteractionRegionCart(CG);
CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);





%% New stff
CG = setupMexInteractionMapping(CG);
[offsets, support, celltypes] = getGridData(CG);
offsets = offsets + 1; 
support = support +1; 

N = G.cells.num;

boundary = zeros(N,1);
     
for (i = 1:length(celltypes))
    if (celltypes(i)==1)
       boundary(support(i))=1; 
    end
end

G.cells.num
G.cells.num / CG.cells.num
sum(boundary)/N

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;


colorMatrix = [ 1 1 1; .55 .55 .55; my_green_2; my_blue_1; my_green_2];


FigHandle = figure('Position', [1200, 200, 500, 500]);
view(50,25);
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,boundary,find(boundary>0.1),'EdgeColor', 'k')
%plotGrid(G,'FaceColor', 'none')
%outlineCoarseGrid(G,pv)
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1
%set(colorbar, 'YTick',1:max(basis+1));


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


