clc; clear all; close all;


nCoarse = 5;
nx = 7*nCoarse; 
ny = nx; 
nBasis = 1; 
nBasis2 = 1; 

G = cartGrid([nx, ny]);
G = computeGeometry(G);

% Create Coarse Grid
pv = partitionUI(G, [nCoarse, nCoarse]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
%CG = storeInteractionRegionCart(CG);
CG = setupMexInteractionMapping(CG);
[offsets, support, celltypes] = getGridData(CG);
offsets = offsets + 1; 
support = support +1; 

N = G.cells.num; 

basis = zeros(1,N);

k = [ 7 8 9 12 13 14 17 18 ]

for (i = 1:length(k))
    for ( j = offsets(k(i)):offsets(k(i)+1)-1 )
        
       if celltypes(j)==2 && basis(support(j)) ~= 3 
           basis(support(j)) = 2;
       elseif basis(support(j)) == 0 && celltypes(j)~=1
           basis(support(j)) = 1; 
       end
       
       if i==5 && (celltypes(j)==2 || celltypes(j)==1 ) 
           basis(support(j)) = 3;
       end
       
       
    end
end

for ( j = offsets(19):offsets(20)-1 )
   if (celltypes(j)==2 && basis(support(j)) == 3)
       support(j)
       basis(support(j)) = 2; 
   end
end


for i = 1:CG.cells.num
   basis(CG.cells.centers(i)) = 5;
end

for i = 1:length(k)
   basis(CG.cells.centers(k(i))) = 3;
end


%basis = zeros(1,N);
for (i = 1:N)
    if (pv(i) > 6 && pv(i)<10 ) || (pv(i) > 11 && pv(i)<15 ) || (pv(i) > 16 && pv(i)<19 )
        %basis(i) = 1; 
    end
end
    







my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;
my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;
my_blue_4 = [34 25 160] ./ 255;
my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;
my_red_3 = [241 36 35] ./ 255;

colorMatrix = [ 1 1 1; .55 .55 .55; my_red_3; my_blue_4; my_green_2];

%colorMatrix = [ 1 1 1; .55 .55 .55];

hold on; 
%plotCellData(G,globalBoundary','EdgeColor', 'y', 'EdgeAlpha',0.1)
plotCellData(G,basis','EdgeColor', 'k')
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'k','linewidth',2)
axis equal tight off; 
%caxis([0 max(pv)]);
colormap(colorMatrix)
colormap(1,1) = 1
%set(colorbar, 'YTick',1:max(basis+1));
zoom(1.1)


%{
q = 9; 
r = 25; 
s = 7
t = -7; 


text (q,r,'a', 'fontsize',30,'color',my_blue_4 ) 
text (q+s,r,'b', 'fontsize',30,'color',my_blue_4 ) 

text (q,r+t,'c', 'fontsize',30,'color',my_blue_4 ) 
text (q+s,r+t,'d', 'fontsize',30,'color',my_blue_4 ) 
text (q+2*s,r+t,'e', 'fontsize',30,'color',my_blue_4 ) 

text (q,r+2*t,'f', 'fontsize',30,'color',my_blue_4 ) 
text (q+s,r+2*t,'g', 'fontsize',30,'color',my_blue_4 ) 
text (q+2*s,r+2*t,'h', 'fontsize',30,'color',my_blue_4 ) 
zoom(1.1)
%}
%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters shareBoundary
%}


%print -dpng temp2


