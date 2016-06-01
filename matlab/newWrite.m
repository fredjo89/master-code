function [] = newWrite(CG, A, varargin)


    opt = struct('writePath', '');
    opt = merge_options(opt, varargin{:});
    
    
    assert(CG.parent.cells.num == size(A, 1));

    
    [offsets, support, types, I_com, fine, coarse, centers] = getGridData(CG);
    
    
    
       
    if ~exist(opt.writePath, 'dir')
        mkdir(opt.writePath)
    end
    inp = fullfile(opt.writePath, 'input');
    if ~exist(inp, 'dir')
    mkdir(inp);
    end
    
    
    [mat, j_index] = compressMatrix(A);
    % Write stuff
    writeGrid(CG, inp, offsets, support, types, centers);
    writeMatrix(inp, j_index', mat');
    writeOperator(inp, I_com);
 
    

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

function writeGrid(CG, fn, offsets, support, types, centers)
    fh = fopen(fullfile(fn, 'info.txt'), 'w');
    fprintf(fh, '%d\r\n%d\r\n', CG.parent.cells.num, CG.cells.num);
    fprintf(fh, '%d ', offsets);
    fprintf(fh, '\n');
    fprintf(fh, '%d ', centers');
    fprintf(fh, '\r\n');
    fclose(fh);

    fh = fopen(fullfile(fn, 'support.txt'), 'w');
    fprintf(fh, '%d ', support');
    fprintf(fh, '\r\n');
    fclose(fh);

    fh = fopen(fullfile(fn, 'types.txt'), 'w');
    fprintf(fh, '%d ', types');
    fprintf(fh, '\r\n');
    fclose(fh);
end

function writeMatrix(fn, jj, mat)
    fh = fopen(fullfile(fn, 'sparsity.txt'), 'w');
    for i = 1:size(mat, 1)
        fprintf(fh, '%d ', jj(i, :));
        fprintf(fh, '\r\n');
    end
    fclose(fh);

    fh = fopen(fullfile(fn, 'matrix.txt'), 'w');
    fprintf(fh, '%d\r\n', size(mat, 2));
    for i = 1:size(mat, 1)
        fprintf(fh, '%1.8f ', mat(i, :));
        fprintf(fh, '\r\n');
    end
    fclose(fh);
end

function writeOperator(fn, I_init)
    fh = fopen(fullfile(fn, 'operator.txt'), 'w');
    fprintf(fh, '%1.8f ', I_init);
    fprintf(fh, '\r\n');
    fclose(fh); 
end