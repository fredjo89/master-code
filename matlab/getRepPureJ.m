function [basis, i] = getRepPureJ(CG, A, iterations, tol, w)




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
    
    tolOne = 10^(-8);
    tolTwo = -1; 
    
    maxIT1 = 10000; 
    maxIT2 = 1; 
    
    
    totalIter = 0; 
     
    j = 0; 
    while j < maxIT1
        j=j+1;
        i = 0;
        prevOne = I; 
        while i < maxIT2
             i = i + 1;
            prevTwo = I; 
            
            
            
            update = M*I;

            update = update.*interactionMap;
            
            
          
            I = I - w*update;

            for k = 1:CG.cells.num
                I(CG.cells.centers(k),k) = 1; 
            end


            diffTwo = max(max(abs(prevTwo - I)));

            if (diffTwo < tolTwo) 
                dispif(false, 'single run converged after %d iterations\n', i);
                break;
            end
            if (i==maxIT2)
                dispif(false, 'single run did not converge after  %d iterations\n', i);
            end
            
        end
        
        totalIter=totalIter+i; 

        
        
        %I = bsxfun(@rdivide, I, sum(I, 2));
        
        diffOne = full(max(max(abs(prevOne - I))));
        if (diffOne < tolOne)
            dispif(true, '\nrepPureJ converged after %d sequences\n', j);
            dispif(true, 'total number of iterations: %d\n', totalIter);
            break;
        end
        
        if (j==maxIT1)
            dispif(true, '\nrepPureJ did not converge after %d sequences\n', j);
            dispif(true, 'total number of iterations: %d\n', totalIter);
        end
        
        
    end
    
    
    
    basis = struct('R', R, 'B', I, 'type', 'rsb');
end


    
