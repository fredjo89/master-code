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

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;

%{
colorMatrix = [ 1 1 1;  my_red_3;  my_green_2];
g = 0.95;
figure(); 
hold on; 
plotCellData(G,unionBoundary,'EdgeColor', 'y', 'EdgeAlpha',0.1)
outlineCoarseGrid(G,pv,'k','lineWidth',3)
plotGrid(G,'FaceColor', 'none')
axis equal tight off; 
colormap(colorMatrix);
%}




for i=1:length(pv)
    if pv(i)~=5
        pv(i) = 0;
    end
end

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters iter0
%}


iterations = 0;



Jbasis = getJBasis(CG,A,iterations, -1, 1);

theBasis = full(Jbasis(:,5));

temp3 = full(unique(Jbasis(:,5)));

for i = 1:length(theBasis)
    for j = 1:length(temp3)
        if theBasis(i)==temp3(j)
            %theBasis(i) = j; 
            break;
        end
    end
end



temp = uniquetol(theBasis);

for i = 1:length(theBasis)
    currentVal = theBasis(i);
    
    
    index = 0; 
    if currentVal ~=0
        index = 1; 
        closest = temp(1); 
        for j = 2:length(temp)
            if abs(closest-currentVal)>abs(temp(j)-currentVal)
                closest = temp(j); 
                index = sqrt(j); 
            end
        end
    end
    
   
   theBasis(i) = index; 

    
end



yolo = pv; 
for (i = 1:length(yolo))
    if theBasis(i)==0
        yolo(i) = 0; 
    else
        yolo(i) = 1; 
    end
end
        

my_green = [0.4 0.4 0.4];


g = 0;
figure(); 
hold on; 
%plotCellData(G,log(theBasis),'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,theBasis,'EdgeColor', 'k', 'EdgeAlpha',0.01)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,support1,'k','linewidth',3.5 )
outlineCoarseGrid(G,yolo,my_green, 'linewidth',3.5 )
outlineCoarseGrid(G,pv,'k', 'linewidth',3.5 )
axis equal tight off; 
colormap(jet(32));

TEMP = jet; 
myColor = [0.65 0.65 0.65]
TEMP(1,:) = [1 1 1];

colormap(TEMP);
%colormap(colorMatrix2);
zoom(1.3);




%print -dpng TEMPFIG1000



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













    
    
    
    