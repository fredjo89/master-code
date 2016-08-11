%% 1D example 
clc; clear all; close all; 




N = 4;                              % Number of cells
n = 2;                                 % Number of cells in each block

M = ceil(N / n);                    % Number of blocks

G = cartGrid([N 1]); 
G = computeGeometry(G);
rock = makeRock(G, 1, 1);

pv = partitionUI(G, [M, 1]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG); 

hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);
gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
initState = initResSol(G,0.0, 1.0); 


full(A)

iterations = 100; 
tol = 10^(-13); 
w = 1;

Jbasis_CONV = getJBasis(CG, A, 100000, tol, 0.95);
Jbasis_CONV = Jbasis_CONV.B; 

Jbasis = getJBasis(CG, A, iterations, -1, 0.95);
Jbasis = Jbasis.B; 


GSbasis = getGSBasis(CG, A, iterations, -1, 1);
GSbasis = GSbasis.B; 





%% Plotting

FigHandle = figure('Position', [1000, 800, 500, 800]);
for i = 1:M
    subplot(M,1,i)
    plot(Jbasis(:,i),'--o'); axis([1 N 0 1]);
end


FigHandle = figure('Position', [2000, 800, 500, 800]);
for i = 1:M
    subplot(M,1,i)
    plot(GSbasis(:,i),'--o'); axis([1 N 0 1]);
end




%bar(Jbasis(:,1))
v = [ 3 4 1 2];

%{
FigHandle = figure('Position', [700, 800, 800, 800]);
for i = 1:M*M
    subplot(M,M,v(i))
    plotCellData(G,full(Jbasis(:,i)));
    outlineCoarseGrid(G,pv,'r')
end

FigHandle = figure('Position', [1500, 800, 800, 800]);
for i = 1:M*M
    subplot(M,M,v(i))
    plotCellData(G,full(GSbasis(:,i)));
    outlineCoarseGrid(G,pv,'r')
end
%}










