% Plots the interaction region of basis function #BNum, and color the
% nodes with different categories
% CATEGORIES: 
% 0: Normal cell, not part of any support boundary. 
% 1: Part of the boundary of the current basis function. 
% 2: Part of the global boundary, but not part of the boundary of the
% current basis function. 
function [] = plotOfCells( CG, BNum )
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
            if celltype(j)==0
                R(i)=10; 
            elseif celltype(j)==1
                R(i) = 30; 
            elseif celltype(j)==2
                R(i)=20; 
            end
            
            j=j+1;
            if j>length(support)
                break;
            end
        end
    end
    plotCellData(G,R); 
    hold on; 
    outlineCoarseGrid(G,CG.partition)


end