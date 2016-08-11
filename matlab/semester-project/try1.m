% Sets up the simplest resivoar model possible and solves it using the multiscale method and the
% incompTPFA function. 

clc; clear all;
close all;

%% SETUP MODEL
[nx ny] = deal(32);
G= cartGrid([nx ny], [500 500]);
G = computeGeometry(G);
rock.perm = ones(G.cells.num,1)*100*milli*darcy;
rock.poro = ones(G.cells.num,1)*.2;


gravity reset off
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);

pv = sum(poreVolume(G,rock)); 

src = addSource([], 1, pv);
src = addSource(src, G.cells.num, -pv);

state = initResSol(G, 0.0, 1.0);

%% CREATE A

hT = computeTrans(G,rock);
hT = hT*1/fluid.properties(state);
hf = G.cells.faces(:,1);
T = 1./accumarray(hf, 1./hT);

nc = G.cells.num; 
i = all(G.faces.neighbors~=0,2);
n1 = G.faces.neighbors(i,1);
n2 = G.faces.neighbors(i,2);
d = accumarray([n1; n2 ], repmat(T(i), [2,1]), [nc, 1]);

I = [ n1; n2 ; (1:nc)']; 
J = [n2 ; n1 ; (1:nc)'];
V = [-T(i); -T(i); d]; clear d; 
A = sparse(double(I), double(J), V,nc,nc);

%% CREATE RHS VECTOR
q = zeros(nc,1);
q(1) = src.rate(1); 
q(nc) = src.rate(2);



%% CREATE CORARSE GRID
pv = partitionUI(G,[4,4]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG);
%plotCellData(G,state.pressure,'EdgeColor', 'y');
%outlineCoarseGrid(G,pv,'k');



%% SOLVE WITH MULTISCALE

basis = getMultiscaleBasis(CG,A);
state_MScale = incompMultiscale(state,CG,hT,fluid,basis,'src',src);
p1 = state_MScale.pressure; 


%% SOLVE WITH REGULAR METHOD
hT = computeTrans(G,rock);
state_2 = incompTPFA(state,G,hT,fluid,'src',src);
p2 = state_2.pressure;

p1=p1*1000;

figure()
plot(A*p1 -q)


sum(abs(A*p2 -q))/length(p2)
figure()
plotCellData(G, state_2.pressure)
plotGrid(G,src.cell)

state_2.pressure = state_2.pressure *10^(-3);
figure()
plot(state_MScale.pressure,'*')
hold on
plot(state_2.pressure,'r*')

