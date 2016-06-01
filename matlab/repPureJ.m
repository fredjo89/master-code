clc; clear all ; close all; 

%% Setting up the model
Temp1 = 4;
Temp2 = 2; 
iterations = 1000; 
tol = 10^(-4);
omega = 2/3;


[nx ny] = deal(Temp1);
theBasisDesider = ceil(Temp1/Temp2);        % Desides how how to devide up the coarse grid.

%% Create Grid
G = cartGrid([nx ny], [2, 2]);
G = computeGeometry(G);

%% Create Coarse Grid
cgxy =ceil(nx/theBasisDesider);
pv = partitionUI(G, [cgxy, cgxy]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);

rock = makeRock(G, 1, 1);

gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
initState = initResSol(G,0.0, 1.0); 

 %% Create A
hT = computeTrans(G,rock);
hT = hT* 1/fluid.properties(initState);
A = getIncomp1PhMatrix(G, hT);


%% Create sourceterms
flowRate = 10^(-16)*sum(poreVolume(G,rock));
src = addSource([], 1,flowRate);
src = addSource(src,G.cells.num ,-flowRate);
for i = 1:nx-1
    src = addSource(src, i*ny+1,flowRate); 
    src = addSource(src, i*ny , -flowRate);
end



Jbasis = getJBasis(CG,A,iterations, tol, omega);
pureBasis = getPureJacobi(CG,A,iterations, tol, omega);
repBasis = getRepPureJ(CG,A,iterations, tol, omega);
repGSBasis = getRepPureGS(CG,A,iterations, tol, omega);



FS_Solution   =  incompTPFA(initState, G, hT, fluid, 'src',src);
J_Solution =  incompMultiscale(initState,CG,hT,fluid,Jbasis,'src',src);
pureSolution =  incompMultiscale(initState,CG,hT,fluid,pureBasis,'src',src);
repSolution = incompMultiscale(initState,CG,hT,fluid,repBasis,'src',src);
repGSSolution = incompMultiscale(initState,CG,hT,fluid,repGSBasis,'src',src);




Jacobi_Diff = abs(FS_Solution.pressure - J_Solution.pressure); 
Jacobi_InfNorm = max(Jacobi_Diff)/max(abs(FS_Solution.pressure))

 
pure_Diff = abs(FS_Solution.pressure - pureSolution.pressure); 
pure_InfNorm = max(pure_Diff)/max(abs(FS_Solution.pressure))

rep_Diff = abs(FS_Solution.pressure - repSolution.pressure); 
rep_InfNorm = max(rep_Diff)/max(abs(FS_Solution.pressure))

repGS_Diff = abs(FS_Solution.pressure - repGSSolution.pressure); 
repGS_InfNorm = max(repGS_Diff)/max(abs(FS_Solution.pressure))


hold on; 
subplot(2,2,1)
plotCellData(G,full( repGSBasis.B(:,3) ));
outlineCoarseGrid(G,pv);

subplot(2,2,2)
plotCellData(G,full( repGSBasis.B(:,4) ));
outlineCoarseGrid(G,pv);

subplot(2,2,3)
plotCellData(G,full( repGSBasis.B(:,1) ));
outlineCoarseGrid(G,pv);

subplot(2,2,4)
plotCellData(G,full( repGSBasis.B(:,2) ));
outlineCoarseGrid(G,pv);



max(max(abs(full(repGSBasis.B-repBasis.B))))


