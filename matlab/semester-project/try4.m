clc; clear all ; close all; 

%% Setting up the model
Temp1 = 4;
Temp2 = 2; 
[nx ny] = deal(Temp1);
theBasisDesider = ceil(Temp1/Temp2);        % Desides how how to devide up the coarse grid.
iterations = 100; 


G = cartGrid([nx ny], [2, 2]);
G = computeGeometry(G);

p = gaussianField(G.cartDims, [.2 .4], [11 3] , 2.5);
K = p.^3*10^(-10)./(0.81*72*(1-p).^2);
rock.poro = p(:); 
rock.perm = K(:); 

%plotCellData(G, rock.perm)
gravity reset off; 
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
pv = sum(poreVolume(G,rock));
src = addSource([], 1,-pv); 
src = addSource(src, G.cells.num, pv); 
initState = initResSol(G,0.0, 1.0); 


 %% Create A
hT = computeTrans(G,rock);
hT = hT* 1/fluid.properties(initState);
A = getIncomp1PhMatrix(G, hT);

full(A)

%% Solve with incompTPFA
%tic
state_fullSol = incompTPFA(initState, G, hT, fluid, 'src',src);
%toc

%% Create Coarse Grid
cgxy =ceil(nx/theBasisDesider);
pv = partitionUI(G, [cgxy, cgxy]);
CG = generateCoarseGrid(G,pv);
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);

%% Get basis
%tic
basis = getMultiscaleBasis(CG,A, 'type', 'rsb', 'iterations',iterations);

state_MScale = incompMultiscale(initState,CG,hT,fluid,basis,'src',src);
%toc



R = ones(size(state_MScale.pressure));
plotCellData(G,state_MScale.pressure);
hold on
outlineCoarseGrid(G,pv)

%{
%% Plotting

figure('Color',[1 1 1]);
R = ones(size(state_MScale.pressure));
plotCellData(G,R)
hold on
outlineCoarseGrid(G,pv)
set(gca,'ycolor',[1,1,1]);
set(gca,'xcolor',[1,1,1]);
set(gca,'ycolor',[1,1,1]);


figure('Color',[1 1 1]);
R = ones(size(state_MScale.pressure))
plotCellData(G,R)
outlineCoarseGrid(G,pv)
set(gca,'ycolor',[1,1,1]);
set(gca,'xcolor',[1,1,1]);
set(gca,'ycolor',[1,1,1]);

B = full(basis.B(:,5)); 

%plotCellData(G,B(:,1));

for k = 1: length(R)
    if B(k,1)~= 0
        B(k,1) = 1; 
    end
end
plotCellData(G,B(:,1));
outlineCoarseGrid(G,pv)

for ii = 1:9
R = full(basis.B(:,ii))
for k = 1: length(R)
    if R(k)~=0
        R(k) =1; 
    end
end
end

subplot(3,3,ii)
plotCellData(G,R)
outlineCoarseGrid(G,pv)
axis equal
%}


%plotCellData(G,state_fullSol.pressure)
%figure()
%plotCellData(G,state_MScale.pressure)
%p_fullsol = state_fullSol.pressure; 
%p_MScale = state_MScale.pressure; 

%subplot(2,2,1)
%{
full(basis.R);
full(basis.B);

plotCellData(G,R)
outlineCoarseGrid(G,pv)

subplot(2,2,2)
plotCellData(G,rock.perm); 
outlineCoarseGrid(G,pv);

subplot(2,2,3)
plotCellData(G,rock.poro); 
outlineCoarseGrid(G,pv);



for k = 1: length(R)
    if R(k)~=0
        R(k) =1; 
    end
end




subplot(2,2,4)
plotCellData(G,R)
outlineCoarseGrid(G,pv)

%figure()
%plot(basis.B(:,1))


%figure()
%plot(p_MS,'r')
%hold on
%plot(p_regular)

norm(A*p_MS - q)
norm(A*p_regular-q)
%}