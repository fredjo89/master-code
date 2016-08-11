% A function that plots the interaction region of 
% basis function #block. 
function [] = plotInteractionRegion( CG, block )
    G = CG.parent; 
    interCells = CG.cells.interaction{block}; 
    ncells = G.cells.num; 
    
    R = zeros(ncells,1); 
    j = 1; 
   
    for i = 1: ncells
        if interCells(j) ==i
            R(i) = 1;
            j=j+1;
            if j>length(interCells)
                break;
            end
        end
    end
    plotCellData(G,R); 
    outlineCoarseGrid(G,CG.partition)
end