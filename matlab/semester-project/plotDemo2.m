% Plotting results

clc; clear all; close all; 
mrstModule add coarsegrid new-multiscale scratchpad libgeometry

%get case
[G, rock, p, testcase] = getCase(); 

if ~isfield(G.cells, 'centroids')
    G = computeGeometry(G);
end

CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG,'edgeBoundaryCenters', false);
%CG = storeInteractionRegionCart(CG);
CG = setupMexInteractionMapping(CG);

%% Read from file
fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);

I2 = cppMultiscaleBasisFromFile_NoMex(fn);


%% Plot of basis functions
mrstModule add mrst-gui
close all;
I2 = full(I2);
%figure; 
hold on; 
R = ones(size(I2,1),1);
basisNumber = 1;
R = I2(:,basisNumber);
%plotCellData(G,R); 
%outlineCoarseGrid(G,p)

initState = initResSol(G,0.0, 1.0); 
hT = computeTrans(G,rock);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
pv = sum(poreVolume(G,rock));
src = addSource([], 1,-pv); 
src = addSource(src, G.cells.num, pv); 
%state_MScale = incompMultiscale(initState,CG,hT,fluid,I2,'src',src);

%figure(); 
plotInteractionRegion(CG,basisNumber);

figure(); 
plotOfCells(CG,basisNumber);

k = 1; 
for i=1:size(I2,2)
    for j = 1:size(I2,1)
        arrayBasis(k) = I2(j,i); 
        k=k+1;
    end
end
figure();
plot(arrayBasis,'x')
norm(arrayBasis)











