clc; clear all; close all;

nx = 4*9; 
ny = 4*9; 
nCoarse = 4;
nBasis = 6; 
nBasis2 = 11; 


G = cartGrid([nx, ny]);
G = computeGeometry(G);


% Create Coarse Grid
pv = partitionUI(G, [nCoarse, nCoarse]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);


currentBasis = CG.cells.interaction{nBasis}; 
currentBasis2 = CG.cells.interaction{nBasis2};
for i=1:nx*ny
    if any(abs(currentBasis-i)<1e-10) && any(abs(currentBasis2-i)<1e-10) && any(abs(pv(i)-11)<1e-10)
        basis(i) = 3; 
    elseif any(abs(currentBasis-i)<1e-10)
        basis(i)=2; 
    elseif any(abs(currentBasis2-i)<1e-10)
        basis(i) = 0; 
    else
        basis(i) = 0; 
    end
    
     if any(abs(CG.cells.centers-i)<1e-10)
        basis(i) = 4; 
    end
    
end

my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_blue_3 = [31 12 146] ./ 255;
my_blue_4 = [34 25 160] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;


%{
%close all; 
globalBoundary = zeros(1,nx*ny);
%find global boundary cells
for i=1:nCoarse*nCoarse
    currentBasis = CG.cells.interaction{i}; 
    currentBoundary = findSupportBoundary(currentBasis,nx,ny);
    for i=1:nx*ny
       if currentBoundary(i)==1
           basis(i)=2;
       end
    end
end
currentBasis = CG.cells.interaction{5}; 
currentBoundary = findSupportBoundary(currentBasis,nx,ny);
for i=1:nx*ny
    if currentBoundary(i)==1
        basis(i)=3;
    end
end


%}



colorMatrix = [ 1 1 1; .35 .35 .35; 0.55 0.55 0.55; my_red_2; my_green_2]


FigHandle = figure('Position', [700, 200, 1300, 1000]);
hold on; 
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotGrid(G,'FaceColor','none','EdgeAlpha',1);
plotCellData(G,basis','EdgeColor', 'k')
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k', 'linewidth',7)
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1
%set(colorbar, 'YTick',1:max(basis+1));

a = 20; 
b = 20; 


text (a,b,'\Omega', 'fontsize',80) 
text (a,b,'j', 'fontsize',40) 
text (a,b+1,'_', 'fontsize',80) 

text (a,b,'\Omega', 'fontsize',80) 
text (a,b,'k', 'fontsize',40) 
text (a,b+1,'_', 'fontsize',80) 




%% ONE
length = 0.09; 
A0 = 0.545; 
A1 = 0.755;
x = [A0 A0];
y = [A1 A1+length];
t = annotation('textarrow',x,y, 'color',my_blue_4)

t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 

xWidth = 0.051;
%1
length = 0.09; 
A0 = 0.545+xWidth; 
A1 = 0.755;
x = [A0 A0];
y = [A1 A1+length];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 
%2
length = 0.09; 
A0 = 0.545+2*xWidth; 
A1 = 0.755;
x = [A0 A0];
y = [A1 A1+length];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 
%3
length = 0.09; 
A0 = 0.545+3*xWidth; 
A1 = 0.755;
x = [A0 A0];
y = [A1 A1+length];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 

%% ONE POINT TWO
length = 0.09; 
A0 = 0.545; 
A1 = 0.465;
x = [A0 A0];
y = [A1+length A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 

%2
A0 = 0.545+xWidth; 
x = [A0 A0];
y = [A1+length A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 
%3
A0 = 0.545+2*xWidth; 
x = [A0 A0];
y = [A1+length A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 
%4
A0 = 0.545+3*xWidth; 
x = [A0 A0];
y = [A1+length A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 6;
t.HeadLength = 25; 
t.HeadWidth = 40; 



%% TWO
length = 0.075; 
A0 = 0.473; 
A1 = 0.755;
x = [A0+length A0];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)

t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40; 

yWidth = 0.067;
%2
length = 0.075; 
A0 = 0.473; 
A1 = 0.755-yWidth;
x = [A0+length A0];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)

t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40; 
%3
length = 0.075; 
A0 = 0.473; 
A1 = 0.755-2*yWidth;
x = [A0+length A0];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)

t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40; 
%4
length = 0.075; 
A0 = 0.473; 
A1 = 0.755-3*yWidth;
x = [A0+length A0];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)

t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40; 




%% TWO POINT 1
length = 0.075; 
A0 = 0.695; 
A1 = 0.755;
x = [A0 A0+length];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40; 
yWidth = 0.067;
% 2
A1 = 0.755-yWidth;
x = [A0 A0+length];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40;
% 3
A1 = 0.755-2*yWidth;
x = [A0 A0+length];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40;
% 4
A1 = 0.755-3*yWidth;
x = [A0 A0+length];
y = [A1 A1];
t = annotation('textarrow',x,y, 'color',my_blue_4)
t.LineWidth = 7;
t.HeadLength = 25; 
t.HeadWidth = 40;




zoom(3);

patch(...
    [2.7 0.24; 2.7 0.36; 2.8 0.3], ...
    [0.05 1.1; -0.05 1.1; 0 1.18], 'k', 'clipping', 'off')












%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters arrowTRY
%}






