function [basis, i] = getPureJacobi(CG, A, iterations, tol, w)

    % Ensuring that the rows of A sum to zero
    A = A - diag(sum(A, 2));
    
    % Creating R the normal way
    R = controlVolumeRestriction(CG.partition);
    
    % Create interaction region 
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);

    G = CG.parent;
    
    % Intiallize basis functions
    I = controlVolumeRestriction(CG.partition)';
    
    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    M = D_inv*A;
        
    i = 0;
    while i < iterations
        last = I; 
        
        update = M*I;
        
        update = update.*interactionMap;
       
        I = I - w*update;
        
        for k = 1:CG.cells.num
            I(CG.cells.centers(k),k) = 1; 
        end
        
        
       % I = bsxfun(@rdivide, I, sum(I, 2));
        
        diff = max(max(abs(last - I)));
        
        if (diff < tol)
            i = i+1; 
            dispif(true, 'pureJacobi converged after %d iterations\n', i);
            break;
        end
        
        i = i + 1;
    end
    
    if (i==iterations)
        dispif(true, 'pureJacobi did not converge after  %d iterations\n', i);
    end
    
    I = bsxfun(@rdivide, I, sum(I, 2));
    
    basis = struct('R', R, 'B', I, 'type', 'rsb');
end


    
