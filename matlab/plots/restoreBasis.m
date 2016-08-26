clc; clear all; close all; 

tic
fileID = fopen('MPI_converged.txt');
C = textscan(fileID,'%f');

MPI_basis = C{1}; 

toc

[G, rock, p, testcase] = makeProblem(); 

toc

G = computeGeometry(G);
T = computeTrans(G, rock);
A = getIncomp1PhMatrix(G, T);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);
CG = setupMexInteractionMapping(CG);

toc



sup = CG.cells.support_mex;
offsets = sup.offsets +1;
celltypes = sup.celltypes;
support = sup.support+1;


I = zeros(G.cells.num, CG.cells.num);
I_0 = controlVolumeRestriction(CG.partition)';

k = 1; 
for i = 1 : CG.cells.num
    for j = offsets(i):offsets(i+1)-1
        if (celltypes(j)~=1)
            I( support(j),i) = MPI_basis(k);
            k = k+1; 
        end
    end
end
toc


a = sum(sum(abs(A*I)));
b = sum(sum(abs(A*I_0)));

a/b;

figure();
plotCellData(G,I(:,1));
outlineCoarseGrid(G,p);

figure();
plotCellData(G,full(I_0(:,1)));
outlineCoarseGrid(G,p);



temp = abs(sum(I')-1);

max(temp)




