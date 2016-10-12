clc; clear all; close all;


nx = 60;
ny = 60;

r=900;

NX = ceil(sqrt(nx*ny/r))

G = cartGrid([nx, ny]);
G = computeGeometry(G);

% Create Coarse Grid
pv = partitionUI(G, [NX, NX]);
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

N
N/CG.cells.num

sum(boundary)/N

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;


colorMatrix = [ 1 1 1; .55 .55 .55; my_green_2; my_blue_1; my_green_2];


hold on; 
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,boundary,find(boundary>0.1),'EdgeColor', 'k')
%plotGrid(G,'FaceColor', 'none')
%outlineCoarseGrid(G,pv,'k','linewidth',3)
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1
%set(colorbar, 'YTick',1:max(basis+1));

%print -dpng 2D_example_100





