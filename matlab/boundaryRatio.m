function [ratio,boundary] = boundaryRatio(CG,G)

CG = setupMexInteractionMapping(CG);

[offsets, support, celltypes] = getGridData(CG);

offsets = offsets +1; 
support = support+1; 

n_sup = length(support); 
n_cells = G.cells.num; 

boundary = zeros(n_cells, 1);

%tic
for (i = 1:n_sup)
   if (celltypes(i)~=0)
       boundary(support(i)) =1;
   end
end
%toc

ratio = sum(boundary)/n_cells;


end