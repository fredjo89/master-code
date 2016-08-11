%% writing testcase to file
 

%% Setting up the model
% get case
[G, rock, p, testcase] = getCase(); 
maxiter = 1; 
tol = 0.001; 

if ~isfield(G.cells, 'centroids')
    G = computeGeometry(G);
end

T = computeTrans(G, rock);
A = getIncomp1PhMatrix(G, T);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG, 'edgeBoundaryCenters', false);
%CG = storeInteractionRegionCart(CG);
CG = setupMexInteractionMapping(CG);

D = diag(A);
n = numel(D);
D_inv = spdiags(1./D, 0, n, n);
full(D_inv*A);



%% Write a testcase to flatfiles on disk

fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)])

cppMultiscaleBasis(CG, A, 'doSolve', false, 'writePath', fn);



