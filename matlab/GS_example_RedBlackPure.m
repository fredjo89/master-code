clc; clear all; close all; 

N = 15;                            % Number of cells
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



support5 = zeros(N*N,1);
redBlack = zeros(N*N,1); 
boundary1 = zeros(N*N,1); 
boundary2 = zeros(N*N,1); 
bothBoundaries = zeros(N*N,1); 
firstUpdate = zeros(N*N,1);
secondUpdate = zeros(N*N,1);
firstSecond = zeros(N*N,1);

for i = offsets(5):offsets(6)-1
    if (celltypes(i)~=1)
        support5(support(i))= 1;
    else
        boundary1(support(i))=1; 
    end
end

for i = 1:N*N
    if redBlack(i)==0
        for j = 1:N*N
            if A(i,j)~=0 && i~=j
                redBlack(j)=1;
            end
        end
    end
    
    if boundary1(i)==1
        for j = 1:N*N
            if A(i,j)~=0 && i~=j && boundary1(j)~=1 && support5(j)~=1
                boundary2(j)=1;
            end
        end
    end
    
end

for i = 1:N*N
    if support5(i)==1
        bothBoundaries(i) = 1; 
    elseif boundary1(i)==1
        bothBoundaries(i) = 2; 
    elseif boundary2(i)==1
        bothBoundaries(i) = 3; 
    end
end

for i = 1:N*N
    if redBlack(i)==1 && (support5(i)==1 || boundary1(i)==1)
        firstUpdate(i) = 1; 
    end
end

for i = 1:N*N
    if redBlack(i)==0 && (support5(i)==1 || boundary1(i)==1 || boundary2(i)==1)
        secondUpdate(i) = 1; 
    end
end


for i = 1:N*N
    if firstUpdate(i)==1
        firstSecond(i) = 1; 
    elseif secondUpdate(i)==1; 
        firstSecond(i)=2;
    end
end

myBlue = [ 0.1 0.1 0.8];
myWhite = [0.9 0.9 0.9];
myGray = [0.55 0.55 0.55];
colorMatrix1 = [1 0 0; 0 0 0  ];

figure(); 
hold on; 
plotCellData(G,redBlack,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,support5,myWhite)
axis equal tight off; 
colormap(colorMatrix1)



colorMatrix1 = [1 1 1 ; myGray;  0.2 0.4 0.2 ; 0.2 0.8  0.2  ];
figure(); 
hold on; 
plotCellData(G,bothBoundaries,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,support5,myWhite)
axis equal tight off; 
colormap(colorMatrix1)



colorMatrix1 = [1 1 1 ; 0 0 0; 1 0 0  ];
figure(); 
hold on; 
plotCellData(G,firstSecond,'EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,support5,myWhite)
axis equal tight off; 
colormap(colorMatrix1)
%}


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters redBlackGS-1
%}




















    
    
    
    