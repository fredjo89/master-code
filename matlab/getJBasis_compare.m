function [error] = getJBasis_compare(CG, A, iterations, w, basisSol)

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
   
    
    basisSol2Norm = sum(basisSol.^2);
    error = zeros(2,iterations+1); 
    
    diff = abs(basisSol - I);
    error(1,1) = max(max(diff));                        % max inf error of the difference matrix
    error(2,1) = max(sqrt(sum(diff.^2)./basisSol2Norm));    % max relative two norm of the difference matrix
    
    for  i =1:iterations
        update = M*I;
        update = update.*interactionMap;
        I = I - w*update;
        I = bsxfun(@rdivide, I, sum(I, 2));
        
        
        diff = abs(basisSol - I);
        error(1,i+1) = max(max(diff));                        % max inf error of the difference matrix
        error(2,i+1) = max(sqrt(sum(diff.^2)./basisSol2Norm));    % max relative two norm of the difference matrix
    end

    
    if (i==iterations)
        dispif(true, 'Jacobi did not converge after  %d iterations\n', i);
    end
    
end


    
