clc; clear all; close all;


[G, rock, p, testcase] = makeProblem(); 

if ~isfield(G.cells, 'centroids')
    G = computeGeometry(G);
end

T = computeTrans(G, rock);
A = getIncomp1PhMatrix(G, T);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);


iterations = 100; 
tol = -1; 
w = 2/3;

[I1, error1] = getJBasis(CG, A, iterations,tol, w);
stop = error1(length(error1)); 
[I2, error2] = getJBasis_lockB(CG, A, iterations,tol, w, stop);





gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

%% Create sourceterms
flowRate = 10^(-16)*sum(poreVolume(G,rock));
src = addSource([], 1,flowRate);
src = addSource(src,G.cells.num ,-flowRate);

for i = 2 : 60
    src = addSource(src, i,flowRate);
    src = addSource(src, G.cells.num -i+1, -flowRate);
end

initState = initResSol(G,0.0, 1.0);

% Fine-scale solver
state_fs = incompTPFA(initState,G,T, fluid, 'src', src); 

% Multiscale solvers
R = controlVolumeRestriction(CG.partition);
basis1 = struct('R', R, 'B', I1, 'type', 'rsb');
basis2 = struct('R', R, 'B', I2, 'type', 'rsb');



state_ms_1 = incompMultiscale(initState,CG,T,fluid,basis1,'src',src);
state_ms_2 = incompMultiscale(initState,CG,T,fluid,basis2,'src',src);

% Compute errorNorms
error = abs(state_fs.pressure - state_ms_1.pressure); 
infNorm = max(error)/max(abs(state_fs.pressure))
twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))


% Compute errorNorms
error = abs(state_fs.pressure - state_ms_2.pressure); 
infNorm = max(error)/max(abs(state_fs.pressure))
twoNorm_error = sqrt(sum(error.^2)/sum(state_fs.pressure.^2))

%{

figure()
plotCellData(G,state_fs.pressure);
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(32)); 
c = colorbar('horiz');  

figure(); 
plotCellData(G,state_ms_1.pressure);
plotGrid(G,src.cell, 'FaceColor', 'w'); 
axis equal tight; colormap(jet(32)); 
c = colorbar('horiz');  

%}




