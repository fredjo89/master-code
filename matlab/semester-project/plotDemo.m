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
CG = storeInteractionRegionCart(CG);
CG = setupMexInteractionMapping(CG);




%% Write a testcase to flatfiles on disk
% We can also write basis functions to disk as a series of text files. This
% is not efficient, but it is useful for setting up and running testcases
% without Matlab.
fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);

disp('Stored files:')
ls(fullfile(fn, 'input'))
%% Pass the filename to the mex function
% The function will read the files and produce output in a folder named
% output. Once it is done, it will be read in and returned as the same type
% of sparse matrix.
I2 = cppMultiscaleBasisFromFile_NoMex(fn);


%% Plot of basis functions
mrstModule add mrst-gui
close all;
I2 = full(I2);
%figure; 
hold on; 
R = ones(size(I2,1),1);
basisNumber = 3;
R = I2(:,basisNumber);
%plotCellData(G,R); 
%outlineCoarseGrid(G,p)


%figure(); 
plotInteractionRegion(CG,basisNumber);

figure(); 
plotOfCells(CG,basisNumber);



