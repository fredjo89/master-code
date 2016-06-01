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
CG = setupMexInteractionMapping(CG);

fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);

newWrite(CG, A, 'writePath', fn);


disp('Stored files:')
ls(fullfile(fn, 'input'))





