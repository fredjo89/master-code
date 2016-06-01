% First script that tests construction of basis functions by using Gauss-
% Seidel iteration. 

clc; clear all ; close all; 

%% Setting up the model
Temp1 = 20;
Temp2 = 2; 
iterations = 1000; 
tol = 10^(-3);


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



%% Get the basis functions
JbasisTwo = getJBasis(CG,A,iterations, tol, 0.95);

GSbasis = getGSBasis(CG,A,iterations, tol, 0.95);

SYMBasis = getSymmetricGSBasis(CG,A,iterations,tol,1.5); 

full(JbasisTwo.B);
full(GSbasis.B);




%% Solving

FSsolution   =  incompTPFA(initState, G, hT, fluid, 'src',src);
JTwoSolution =  incompMultiscale(initState,CG,hT,fluid,JbasisTwo,'src',src);
GSsolution   =  incompMultiscale(initState,CG,hT,fluid,GSbasis,'src',src);
SYMsolution   =  incompMultiscale(initState,CG,hT,fluid,SYMBasis,'src',src);

%% Computing discrepancy norms


JacobiTwo_Diff = abs(FSsolution.pressure - JTwoSolution.pressure); 
JacobiTwo_InfNorm = max(JacobiTwo_Diff)/max(abs(FSsolution.pressure));
JacobiTwo_twoNorm = sqrt(sum(JacobiTwo_Diff.^2)/sum(FSsolution.pressure.^2));

GS_Diff = abs(FSsolution.pressure - GSsolution.pressure); 
GS_InfNorm = max(GS_Diff)/max(abs(FSsolution.pressure));
GS_twoNorm = sqrt(sum(GS_Diff.^2)/sum(FSsolution.pressure.^2));

SYM_Diff = abs(FSsolution.pressure - SYMsolution.pressure); 
SYM_InfNorm = max(SYM_Diff)/max(abs(FSsolution.pressure));
SYM_twoNorm = sqrt(sum(SYM_Diff.^2)/sum(FSsolution.pressure.^2));


JacobiTwo_InfNorm
GS_InfNorm
SYM_InfNorm



hold on; 
subplot(2,2,1)
plotCellData(G,full( SYMBasis.B(:,1) ));
outlineCoarseGrid(G,pv);

subplot(2,2,2)
plotCellData(G,full( JbasisTwo.B(:,1) ));
outlineCoarseGrid(G,pv);

subplot(2,2,3)
plotCellData(G,full( SYMBasis.B(:,4) ));
outlineCoarseGrid(G,pv);

subplot(2,2,4)
plotCellData(G,full( JbasisTwo.B(:,4) ));
outlineCoarseGrid(G,pv);


sum(full(SYMBasis.B-JbasisTwo.B))







%% Plotting

%{

hold on; 
subplot(2,2,1)
plotCellData(G,FSsolution.pressure);
outlineCoarseGrid(G,pv);


subplot(2,2,3)
plotCellData(G,JTwoSolution.pressure);
outlineCoarseGrid(G,pv);

subplot(2,2,4)
plotCellData(G,GSsolution.pressure);
outlineCoarseGrid(G,pv);

figure();
hold on; 
subplot(2,2,1)
plotCellData(G,FSsolution.pressure);
outlineCoarseGrid(G,pv);


subplot(2,2,3)
plotCellData(G,JTwoSolution.pressure-FSsolution.pressure);
outlineCoarseGrid(G,pv);

subplot(2,2,4)
plotCellData(G,GSsolution.pressure-FSsolution.pressure);
outlineCoarseGrid(G,pv);


figure(); 
subplot(1,2,1)
plot(JTwoSolution.pressure-FSsolution.pressure,'*')
subplot(1,2,2)
plot(GSsolution.pressure-FSsolution.pressure,'*')
%}

