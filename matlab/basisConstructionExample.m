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

for i=1:length(pv)
    if pv(i)~=5
        pv(i) = 0;
    end
end





for i = 1:SIZE
   if celltypes(i)==1
       unionBoundary(support(i))=1; 
   end
   if i>offsets(5) &&i<offsets(6)
       if (celltypes(i)~=1)
        support1(support(i)) = 1; 
       else
          % support1(support(i)) = 2;
       end
   end
end

for i = 1:CG.cells.num
    unionBoundary(CG.cells.centers(i)) = 3; 
end


iterations = 1000;


Jbasis = getJBasis(CG,A,iterations, -1, 1);

theBasis = full(Jbasis.B(:,5));

temp3 = full(unique(Jbasis.B(:,5)));

for i = 1:length(theBasis)
    for j = 1:length(temp3)
        if theBasis(i)==temp3(j)
           % theBasis(i) = j; 
            break;
        end
    end
end



g = 0.95;
figure(); 
hold on; 
plotCellData(G,log(theBasis),'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,support1,[g g g])
outlineCoarseGrid(G,pv,[g g g])
axis equal tight off; 
colormap(jet(100000));
%colormap(colorMatrix2);


set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters iter1000

%{
parula	
jet	
hsv	
hot	
cool	
spring	
summer	
autumn	
winter	
gray	
bone	
copper	
pink	
lines	
colorcube	
prism	
flag	
white
%}





%{
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
print -dpdf -painters iter0
%}













    
    
    
    