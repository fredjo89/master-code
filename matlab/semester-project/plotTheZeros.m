% Plots the zero valued cells of bais #block
function [] = plotTheZeros( CG, BNum )
    G = CG.parent; 
    ncells = G.cells.num; 
    
    c_total=CG.cells.support_mex.celltypes; 
    s_total = CG.cells.support_mex.support;
    
    support = s_total(CG.cells.support_mex.offsets(BNum)+1:CG.cells.support_mex.offsets(BNum+1))+1;
    celltype = c_total(CG.cells.support_mex.offsets(BNum)+1:CG.cells.support_mex.offsets(BNum+1));
    
    R = zeros(ncells,1); 
    j = 1; 
    for i = 1: ncells
        currentCell = support(j);
        if currentCell == i
            if celltype(j)==1
                R(i)=1; 
            end
            j=j+1;
            if j>length(support)
                break;
            end
        end
    end
    plotCellData(G,R); 
    %outlineCoarseGrid(G,CG.partition)
    
    
    
end