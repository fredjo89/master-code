clc; clear all; close all; 

N = 3*9;                            % Number of cells
n = N/3;                             % Number of cells in each block
M = ceil(N / n);                    % Number of blocks
G = cartGrid([N N]); 
G = computeGeometry(G);
rock = makeRock(G, 1, 1);
pv = partitionUI(G, [M, M]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
%CG = storeInteractionRegionCart(CG);
CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);
hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);




% Create interaction region 
lens = cellfun(@numel, CG.cells.interaction);
blocks = rldecode((1:CG.cells.num)', lens);
interaction = vertcat(CG.cells.interaction{:});
interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);
I = controlVolumeRestriction(CG.partition)';
    
[offsets, support, celltypes] = getGridData(CG);
offsets = offsets + 1; 
support = support + 1; 
SIZE = length(support); 


unionBoundary = zeros(N*N,1);
support1 = zeros(N*N,1);
support9 = zeros(N*N,1);

for i = 1:SIZE
   if celltypes(i)==1
       unionBoundary(support(i))=1; 
   end
   if i<offsets(2)
       if (celltypes(i)~=1)
        support1(support(i)) = 1; 
       else
           support1(support(i)) = 2;
       end
   elseif offsets(length(offsets)-1)<=i
       if (celltypes(i)~=1)
        support1(support(i)) = 1; 
       else
           support1(support(i)) = 2;
       end
   end
end

for i = 1:CG.cells.num
    support1(CG.cells.centers(i)) = 3; 
    unionBoundary(CG.cells.centers(i)) = 3; 
end


iterations = 1;

Jbasis = getJBasis(CG,A,iterations, -1, 0.95);

GSbasis = getGSBasis(CG,A,iterations, -1, 0.95);


%colorMatrix = [ 1 1 1; .55 .55 .55; 1 .3 0; 0 .3 1; 0 0.75 0];

colorMatrix1 = [ 1 1 1; 1 0.2 0; 0 0.55 0 ];
figure(); 
hold on; 
plotCellData(G,unionBoundary,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
colormap(colorMatrix1)
%set(colorbar, 'YTick',1:max(basis+1));


colorMatrix2 = [ 1 1 1; 0.55 0.55 0.55; 1 0.2 0; 0 0.55 0];
figure(); 
hold on; 
plotCellData(G,support1,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
colormap(colorMatrix2)



figure(); 
hold on; 
plotCellData(G,log(full(GSbasis.B(:,1))),'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
colormap(jet(1000));


figure(); 
hold on; 
plotCellData(G,log(full(GSbasis.B(:,9))),'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
colormap(jet(1000));

%{
temp = [1:N*N]';

figure(); 
hold on; 
plotCellData(G,temp,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
colormap(jet(32*2*2*2*2*2*2));
colorbar('horiz')
%}

%{
figure(); 
hold on; 
plotCellData(G,log(full(Jbasis.B(:,1))),'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 

figure(); 
hold on; 
plotCellData(G,full(Jbasis.B(:,9)),'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k')
axis equal tight off; 
%}


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters numbering
%}













    
    
    
    