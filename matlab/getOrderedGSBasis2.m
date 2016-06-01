function [basis, i] = getOrderedGSBasis2(CG, A, iterations, tol, w)

    A = A - diag(sum(A, 2));
    
    CG = setupMexInteractionMapping(CG);
    [offsets, support, types, I_com, fine, coarse, centers] = getGridData(CG);
     
    [mat, j_index] = compressMatrix(A);
    
    Orders = getOrders(CG,A);
    
    
    
    
    
    
    
    
    
    
    
    
    R = controlVolumeRestriction(CG.partition);
    
    I = controlVolumeRestriction(CG.partition)';
    i = 1; 
    basis = struct('R', R, 'B', I, 'type', 'rsb');
end





function [offsets, support, celltypes, I, fine, coarse, centers] = getGridData(CG)
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
    
    
    centers = CG.cells.centers -1;
    
end

function [mat, jj] = compressMatrix(A)
    n = size(A, 1);
    
    D = spdiags(1./diag(A), 0, n, n);
    [i_ix, j_ix, av] = find((D*A)');
    intx = i_ix ~= j_ix;
    i_ix = i_ix(intx);
    j_ix = j_ix(intx);
    av = av(intx);

    j2 = j_ix;

    pos = 1;
    t = tic();
    for i = 1:n
        v = j_ix(pos);
        ctr = 1;
        while v == i
            j2(pos) = ctr;
            ctr = ctr + 1;
            pos = pos + 1;
            if pos > numel(j_ix)
                break
            end
            v = j_ix(pos);
        end
    end
    
    mat = full(sparse(j_ix, j2, av))';
    jj = toIntegerIndex(full(sparse(j_ix, j2, i_ix)))';

    toc(t);
end

function v = toIntegerIndex(v)
    v = toIntegerValue(v - 1);
end

function v = toIntegerValue(v)
    v = int32(v);
end




