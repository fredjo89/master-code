clc; clear all; close all; 

% Number of cells
nx = 60;
ny = 30; 
nz = 10; 
% Number of coarse blocks
NX = 5;
NY = 3;
NZ = 2;

G = processGRDECL(simpleGrdecl([nx ny nz], 0.0));
G = computeGeometry(G);
%rock.perm = logNormLayers(G.cartDims, [100 400 50 350], 'indices', [1 2 5 7 11]);

%rock.perm = ones(G.cells.num,1);

makeRock = 0; 

if makeRock==0
    RR = 6; 
    rock.perm = logNormLayers(G.cartDims, [1000 ], 'indices', [1 11], 'std', 4.5, 'sigma',10, 'a',0.6, 'sz', [15+RR,7+RR,5+RR]);
    fid=fopen('MyRock.txt','w');
    fprintf(fid, '%f \n', rock.perm);
    fclose(fid); 
else
    fileID = fopen('MyRock.txt');
    C = textscan(fileID,'%f');
    rock.perm = C{1}; 
end
%rock.perm = convertFrom(rock.perm,milli*darcy);



gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

% Boundary conditions 
bc = pside([], G, 'left', 1);
bc = pside(bc, G, 'right', 0);

initState = initResSol(G,0.0, 1.0);
hT = computeTrans(G,rock); 

% Fine-scale solver
state_fs = incompTPFA(initState,G,hT, fluid,  'bc', bc); 

% MsRSB solver
A = getIncomp1PhMatrix(G, hT);
%p = partitionUI(G,[NX NY NZ]); 
p = partitionUniformPadded(G, [NX NY NZ]);
G_fault = makeInternalBoundary(G, find(G.faces.tag > 0))
p = processPartition(G_fault, p);

CG = generateCoarseGrid(G,p); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegion(CG); 
CG = setupMexInteractionMapping(CG);

x = [ 0 2 4 7 10 15 20 25 30 40 50 60 70 80 90 100]; 
%x = 0;
for i = 1: length(x)
basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',x(i));

basisRB = basis; 

basisRB.B = getRedBlackBasis(CG, A, x(i),-1, 1);



state_ms = incompMultiscale(initState,CG,hT,fluid,basis, 'bc',bc);

state_ms_RB = incompMultiscale(initState,CG,hT,fluid,basisRB, 'bc',bc);


% Calculating error norms
error = abs(state_fs.pressure - state_ms.pressure); 
infNorm(i) = max(error)/max(abs(state_fs.pressure));
twoNorm_error(i) = sqrt(sum(error.^2)/sum(state_fs.pressure.^2));



error = abs(state_fs.pressure - state_ms_RB.pressure); 
infNorm_2(i) = max(error)/max(abs(state_fs.pressure));
twoNorm_error_2(i) = sqrt(sum(error.^2)/sum(state_fs.pressure.^2));


twoNorm_error
twoNorm_error_2

end



figure
plot(x,twoNorm_error)
hold on; 
plot(x,twoNorm_error_2)

figure
plot(x,infNorm)
hold on; 
plot(x,infNorm_2)



%{
rock.perm = ones(G.cells.num,1);
hT = computeTrans(G,rock); 
state_fs_2 = incompTPFA(initState,G,hT, fluid,  'bc', bc); 
error = abs(state_fs.pressure - state_fs_2.pressure); 
max(error)/max(abs(state_fs.pressure))
sqrt(sum(error.^2)/sum(state_fs.pressure.^2))
%}








