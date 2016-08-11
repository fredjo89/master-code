% A new write function that removes type 1 cells. 

function [] = newWrite2(CG, A, varargin)


    opt = struct('writePath', '');
    opt = merge_options(opt, varargin{:});
    
    
    assert(CG.parent.cells.num == size(A, 1));

    
    [offsets, support, types, I_com, fine, coarse] = getGridData(CG);
    
    n_ones = sum(types == 1);
    n_total = offsets(length(offsets)); 
    n_new = n_total - n_ones; 
    
    offsets_new = zeros(length(offsets),1);
    support_new = zeros(n_new,1);
    types_new = zeros(n_new,1);
    I_com_new = zeros(n_new,1);
    
    k = 1; 
    supportNumber = 1; 
    for i = 1 : n_total
        if types(i)~=1
            
            if i>offsets(supportNumber+1)
                supportNumber=supportNumber+1;
            end
            
            offsets_new(supportNumber+1) = offsets_new(supportNumber+1)+1; 
            support_new(k) = support(i); 
            types_new(k) = types(i); 
            I_com_new(k) = I_com(i); 
            k=k+1; 
        end
    end
    
    for (i=3:length(offsets))
        offsets_new(i)=offsets_new(i)+offsets_new(i-1);
    end
        
    offsets = offsets_new; 
    support = support_new; 
    types   = types_new; 
    I_com   = I_com_new; 

       
    if ~exist(opt.writePath, 'dir')
        mkdir(opt.writePath)
    end
    inp = fullfile(opt.writePath, 'input');
    if ~exist(inp, 'dir')
    mkdir(inp);
    end
    
    
    [mat, j_index] = compressMatrix(A);
    % Write stuff
    writeGrid(CG, inp, offsets, support, types);
    writeMatrix(inp, j_index', mat');
    writeOperator(inp, I_com);
 
    

end

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

function writeGrid(CG, fn, offsets, support, types)
    fh = fopen(fullfile(fn, 'info.txt'), 'w');
    fprintf(fh, '%d\r\n%d\r\n', CG.parent.cells.num, CG.cells.num);
    fprintf(fh, '%d ', offsets);
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