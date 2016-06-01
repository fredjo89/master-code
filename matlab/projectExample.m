%% Homo case
clc; clear all; close all; 

nBlocksX = 5;
nBlocksY = nBlocksX; 
nCellsX = 75; 
nCellsY = nCellsX; 
lengthX = 500; 
lengthY = lengthX;
midWidth = 15;
channelWidth = 7; 

theBasis = 12; 

outerPoro = [ 0.3 0.5];
channelPoro = [ 0.6 0.65];
innerPoro = [0.1 0.2];
iterations = 600; 

% makeRock = 1 if we want to create new rock fields. Otherwise, it reads
makeRock = 1; 

G = cartGrid([nCellsX, nCellsY], [lengthX, lengthY]);
G = computeGeometry(G);

if makeRock==1
    rock = getRock( G, nCellsX, nCellsY,midWidth, channelWidth, outerPoro,... 
                    innerPoro, channelPoro );
    fid=fopen('MyRock.txt','w');
    fprintf(fid, '%f \n', rock.poro);
    fclose(fid); 
else
    fileID = fopen('MyRock.txt');
    C = textscan(fileID,'%f');
    rock.poro = C{1}; 
    rock.perm = rock.poro.^3*10^(-10)./(0.81*72*(1-rock.poro).^2);
end
            
gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);


%% Create sourceterms
flowRate = 10^(-16)*sum(poreVolume(G,rock));
src = addSource([], 1,flowRate);
src = addSource(src,G.cells.num ,-flowRate);
for i = 1:nCellsX-1
    src = addSource(src, i*nCellsY+1,flowRate); 
    src = addSource(src, i*nCellsY , -flowRate);
end


initState = initResSol(G,0.0, 1.0);

hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);

% Fine-scale solver
state_fs = incompTPFA(initState,G,hT, fluid, 'src', src); 

% Creating CG
pv = partitionUI(G,[nBlocksX,nBlocksY]); 
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG); 

% MsRSB solver
basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',iterations);
state_ms = incompMultiscale(initState,CG,hT,fluid,basis,'src',src);


% Compute errorNorms
error = abs(state_fs.pressure - state_ms.pressure); 
infNorm = max(error)/max(abs(state_fs.pressure))
twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))

%% Plotting results 
lWidth = 1; 
fSize = 50; 
axisSize = 35; 


figure(); 
plotCellData(G,state_fs.pressure);
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(32)); 
c = colorbar('horiz');  
set(gca, 'FontSize',axisSize);

figure(); 
plotCellData(G,state_ms.pressure);
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(32)); 
c=colorbar('horiz'); 
set(gca, 'FontSize',axisSize);


figure(); 
plotCellData(G,rock.poro);
outlineCoarseGrid(G,pv,'r')
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(1000*128)); 
c=colorbar('horiz'); 
set(gca, 'FontSize',axisSize);

figure(); 
plotGrid(G,'FaceColor', 'none')
outlineCoarseGrid(G,pv,'r')
axis tight off; 
%caxis([0 max(pv)]);
colormap(1,1) = 1

%{
figure(); 
oneBasis = full(basis.B(:,theBasis)); 
plotCellData(G,oneBasis);
outlineCoarseGrid(G,pv)
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(128*2)); 
c=colorbar('horiz'); 
set(gca, 'FontSize',axisSize);
%}

%{
plot(rock.poro,'x')
figure(); 
plot(rock.perm,'x')
figure(); 
plot(log(rock.perm),'x')
%}

%plot(state_fs.pressure)


%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters epsFig
%}

