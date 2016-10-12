clc; clear all; close all; 



%% Create triangluar plot
N=7; M=5; [x,y] = ndgrid(0:N, 0:M); 
x(2:N,2:M) = x(2:N,2:M)+0.3*randn(N-1, M-1);
y(2:N,2:M) = x(2:N,2:M)+0.3*randn(N-1, M-1) ;
aG = pebi(triangleGrid([x(:) y(:)]));
G = makeLayeredGrid(aG,3); 
plotGrid(G,'FaceColor', [.8 .8 .8]); view(-40,60); axis tight off

clear all; close all; 

N = 21; 
M = 31; 
[x,y] = meshgrid(0:M-1, 0:N-1);


x(1:N, 2:M-1) = x(1:N,2:M-1)+0.2*randn(N, M-2);
y(2:N-1, 1:M) = y(2:N-1,1:M)+0.2*randn(N-2, M);

t = delaunay(x(:),y(:)); 
G = triangleGrid([x(:) y(:)],t);
G = computeGeometry(G);

for i=1:G.cells.num
    x = G.cells.centroids(i,1); 
    y = G.cells.centroids(i,2); 
    if x<8
        if y<7
            pv(i,1) =1; 
        elseif y<13
            pv(i,1) =2; 
        else
            pv(i,1)=3; 
        end
    elseif x<16
        if y<7
            pv(i,1) =4; 
        elseif y<13
            pv(i,1) =5; 
        else
            pv(i,1)=6; 
        end
    elseif x<24
        if y<7
            pv(i,1) =7; 
        elseif y<13
            pv(i,1) =8; 
        else
            pv(i,1)=9; 
        end
    else 
        if y<7
            pv(i,1) =1; 
        elseif y<13
            pv(i,1) =2; 
        else
            pv(i,1)=3; 
        end
    end
end

%% Create Coarse Grid
%pv = partitionUI(G, [2, 2]);
%CG = generateCoarseGrid(G,pv);
%CG = coarsenGeometry(CG); 
%CG = storeInteractionRegionCart(CG);


r = 1; 
FigHandle = figure('Position', [500, 500, 1000*r, 500*r ]);
plotCellData(G,pv,'EdgeColor', 'y')
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k','linewidth',3)
axis tight off; 
caxis([.5 max(pv)+.5]);
colormap(lines(max(pv)));


%set(colorbar, 'YTick',1:max(pv));












