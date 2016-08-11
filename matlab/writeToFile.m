clc; clear all; close all; 


disp('makeProblem:')
tic
[G, rock, p, testcase] = makeProblem(); 
toc

disp('computeGeometry:')
tic
if ~isfield(G.cells, 'centroids')
    G = computeGeometry(G);
end
toc

disp('computeTrans:')
tic 
T = computeTrans(G, rock);
toc

disp('getIncomp1PhMatrix:')
tic
A = getIncomp1PhMatrix(G, T);
toc

disp('generateCoarseGrid:')
tic
CG = generateCoarseGrid(G, p);
toc

disp('coarsenGeometry:')
tic 
CG = coarsenGeometry(CG);
toc

disp('storeInteractionRegionCart:')
tic 
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
toc

disp('setupMexInteractionMapping:')
tic
CG = setupMexInteractionMapping(CG);
toc

fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);
%fn = fullfile('/global/work/fredjoha/mrst_data', ['basis_', lower(testcase)]);



%cppMultiscaleBasis(CG, A, 'doSolve', false, 'writePath', fn);
newWrite2(CG, A, 'writePath', fn);


disp('Stored files:')
ls(fullfile(fn, 'input'))

%full(getJBasis(CG, A, 100,0.00001, 2/3))


temp = zeros(G.cells.num,1);

for i = 1:CG.cells.num
    temp(CG.cells.centers(i)) =1;
end


    %plotCellData(G,temp);
    %outlineCoarseGrid(G,p,'r')


