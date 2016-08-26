clc; clear all; close all; 

fileID = fopen('MPI_converged.txt');
C = textscan(fileID,'%f');

MPI_basis = C{1}; 

[G, rock, p, testcase] = makeProblem(); 

G = computeGeometry(G);
T = computeTrans(G, rock);
A = getIncomp1PhMatrix(G, T);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);


Jbasis = getJBasis(CG, A, 20000, -1, 2/3);
Jbasis = full(Jbasis);


sup = CG.cells.support_mex;
offsets = sup.offsets+1;
celltypes = sup.celltypes;
support = sup.support+1;




MPIbasis = zeros(G.cells.num, CG.cells.num);
I_0 = controlVolumeRestriction(CG.partition)';

k = 1; 
for i = 1 : CG.cells.num
    for j = offsets(i):offsets(i+1)-1
        if (celltypes(j)~=1)
            MPIbasis( support(j),i) = MPI_basis(k);
            k = k+1; 
        end
    end
end
toc


max(max(abs(Jbasis-MPIbasis)))

MPIbasis = bsxfun(@rdivide, MPIbasis, sum(MPIbasis, 2));


max(max(abs(Jbasis-MPIbasis)))





a = sum(sum(abs(A*MPIbasis)));
b = sum(sum(abs(A*Jbasis)));
c = sum(sum(abs(A*I_0)));


a/c;

b/c;

%{
figure();
plotCellData(G,full(I_0(:,1)));
outlineCoarseGrid(G,p);


figure();
plotCellData(G,full(Jbasis(:,27)));
outlineCoarseGrid(G,p);

figure();
plotCellData(G,full(MPIbasis(:,27)));
outlineCoarseGrid(G,p);

a = sum(sum(abs(A*Jbasis)));
b = sum(sum(abs(A*I_0)));

a/b


a = sum(sum(abs(A*MPIbasis)));
b = sum(sum(abs(A*I_0)));

a/b

%}






















