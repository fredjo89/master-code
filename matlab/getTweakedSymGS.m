function [basis, i] = getTweakedSymGS(CG, A, iterations, tol, w)


    iterations = floor(iterations);

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
    
    RHS = sparse(G.cells.num, CG.cells.num);
    
    
    
    for k = 1:CG.cells.num
        temp = zeros(1,G.cells.num); 
        temp(CG.cells.centers(k)) = 1; 
        A(CG.cells.centers(k),:) =temp; 
        
        RHS(CG.cells.centers(k),k)=1;
    end
    
    
    F = - (triu(A) - spdiags(D, 0, n, n));
    E = - (tril(A) - spdiags(D, 0, n, n)) ;
    
  
    FORWARD = A+F;
    BACKWARD = A+E;
    
    tol = tol; 

    
    maxIT = iterations;  
    
    totalIter = 0; 
   
        
    j = 0; 
    while j < maxIT
        j=j+1;
        prev = I; 
        
        
        if rem(j,2)==1
            % forward step
            update = (FORWARD)\(F*I+RHS); 
            update = update.*interactionMap;
            update = (1-w)*I + w*update;
            I = update;
            I = bsxfun(@rdivide, I, sum(I, 2));
        else

            %backward step
            update = BACKWARD\(E*I+RHS); 
            update = update.*interactionMap;
            update = (1-w)*I + w*update;
            I = update;
            I = bsxfun(@rdivide, I, sum(I, 2));
        end
        
        diff = full(max(max(abs(prev - I))));
        if (diff < tol)
            dispif(true, '\nSymmetricGS converged after %d sequences\n', j);
            dispif(true, 'total number of iterations: %d\n', j);
            break;
        end
        
        if (j==maxIT)
            dispif(true, '\nSymmetricGS did not converge after %d sequences\n', j);
            dispif(true, 'total number of iterations: %d\n', j);
        end
        
        
    end
    
    i = j; 
    

    
    basis = struct('R', R, 'B', I, 'type', 'rsb');
end


    
