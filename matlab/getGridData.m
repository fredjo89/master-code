function [offsets, support, celltypes, I, fine, coarse] = getGridData(CG)
    sup = CG.cells.support_mex;
    offsets = sup.offsets;
    celltypes = sup.celltypes;
    support = sup.support;
    
    % Convert to 1-indexing
    ofs = double(offsets) + 1;
    I = zeros(size(support));
    for i = 1:CG.cells.num
        subs = ofs(i):(ofs(i+1)-1);
        I(subs) = CG.partition(support(subs, 1) + 1) == i;
    end
    
    fine = vertcat(sup.sorted_cell{:});
    coarse = rldecode((1:CG.cells.num)', cellfun(@numel, sup.sorted_cell));

end