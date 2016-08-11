clc; clear all; close all; 

fileID = fopen('MPI_basis.txt');
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


Jbasis = getJBasis(CG, A, 100, -1, 2/3);


sup = CG.cells.support_mex;
offsets = sup.offsets;
celltypes = sup.celltypes;
support = sup.support;



Jbasis_comprs = zeros(length(MPI_basis),1);
k = 1; 
for i = 1:CG.cells.num    
    for j = offsets(i)+1:offsets(i+1)
       if celltypes(j)~=1
           Jbasis_comprs(k) = Jbasis(support(j)+1,i);
           k = k+1; 
       end
    end
    
end

sum(abs(Jbasis_comprs-MPI_basis))/length(MPI_basis)


plotCellData(G,full(Jbasis(:,2)));
outlineCoarseGrid(G,p,'r')



