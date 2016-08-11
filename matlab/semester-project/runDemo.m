% computing bassis function, no plotting. 
 

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
CG = storeInteractionRegionCart(CG);
CG = setupMexInteractionMapping(CG);

D = diag(A);
n = numel(D);
D_inv = spdiags(1./D, 0, n, n);
full(D_inv*A);

%% Generate basis functions via regular mex interface
% Typical usage. We send in the coarsegrid and the mex layer passes values
% onto the C++ code.
I = cppMultiscaleBasis(CG, A, 'verbose', true, 'omega', 2/3, 'maxiter', maxiter, 'tolerance', tol);



%% Write a testcase to flatfiles on disk
% We can also write basis functions to disk as a series of text files. This
% is not efficient, but it is useful for setting up and running testcases
% without Matlab.

fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);
rmdir(fn,'s');
cppMultiscaleBasis(CG, A, 'doSolve', false, 'writePath', fn);

%disp('Stored files:')
%ls(fullfile(fn, 'input'))
%% Pass the filename to the mex function
% The function will read the files and produce output in a folder named
% output. Once it is done, it will be read in and returned as the same type
% of sparse matrix.
I2 = cppMultiscaleBasisFromFile(fn);



