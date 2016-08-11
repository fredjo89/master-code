function [basis, i] = getOrderedGSBasis(CG, A, iterations, tol, w)
 % Setup
    Orders = getOrders(CG,A);
    A = A - diag(sum(A, 2));
    R = controlVolumeRestriction(CG.partition);
    
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);

    I = controlVolumeRestriction(CG.partition)';
    
    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    
   M = D_inv*A;
    
    C = CG.cells.centers;
    
    i = 0;
    while i < iterations
        update = I; 
        % get update
        for (n_b = 1:CG.cells.num)
            for (n_c_index = 2:n)
                r = Orders(n_c_index, n_b);
                update(r,n_b) = 0; 
                for (j = 1:n)
                    if (j~=r)
                        update(r,n_b) = update(r,n_b) - M(r,j)*update(j,n_b);
                    end
                end
            end
        end
       
        update = update.*interactionMap;
        update = bsxfun(@rdivide, update, sum(update, 2));
        
        diff = max(max(abs(update-I)));
        I = update; 
        
        if (diff < tol)
            i = i + 1;
            dispif(true, 'Orered GS converged after %d iterations\n', i);
            break;
        end
        i = i + 1;
    end
    
     
    

    dispif(i==iterations, 'Ordered GS DID NOT converged after %d iterations\n', i);
    basis = struct('R', R, 'B', I, 'type', 'rsb');
end



