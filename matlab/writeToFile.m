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


fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);
%fn = fullfile('/global/work/fredjoha/mrst_data', ['basis_', lower(testcase)]);



%cppMultiscaleBasis(CG, A, 'doSolve', false, 'writePath', fn);
newWrite2(CG, A, 'writePath', fn);


disp('Stored files:')
ls(fullfile(fn, 'input'))

%full(getJBasis(CG, A, 100,0.00001, 2/3))

%{
temp = zeros(G.cells.num,1);

for i = 1:CG.cells.num
    temp(CG.cells.centers(i)) =1;
end
plotCellData(G,temp);
outlineCoarseGrid(G,p,'k')
%}

