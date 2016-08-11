function [basis, i] = getGSBasis(CG, A, iterations, tol, w)
    A = A - diag(sum(A, 2));
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
    
    F = - (triu(A) - spdiags(D, 0, n, n));
    E = - (tril(A) - spdiags(D, 0, n, n)) ;
    
  
    FORWARD = A+F;
    BACKWARD = A+E;

    i = 0;
    while i < iterations
        
        update = (FORWARD)\(F*I); 
        %update = BACKWARD\(E*I); 
        
        % forward step
        %update = (FORWARD)\(F*I); 
        %update = update.*interactionMap;
        %update = bsxfun(@rdivide, update, sum(update, 2));
      
        %update = (1-w)*I + w*update;
        
        %backward step
        %update = BACKWARD\(E*update); 
        %update = update.*interactionMap;
        %update = bsxfun(@rdivide, update, sum(update, 2));
  
        %update = (1-w)*I + w*update;
        
        
        
        
        diff = max(max(abs(update-I)));
        
        I = update; 
        
        if (diff < tol)
            i = i + 1;
            dispif(true, 'Gauss-Seidiel converged after %d iterations\n', i);
            break;
        end
        
        i = i + 1;
    end
    
    if (i==iterations)
        dispif(true, 'Gauss-Seidiel DID NOT converged after %d iterations\n', i);
    end
    
    basis = struct('R', R, 'B', I, 'type', 'rsb');
end



