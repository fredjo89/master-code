function [basis, i] = getRepPureGS(CG, A, iterations, tol, w)


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
  
    FORWARD = A+F;

    tolOne = tol; 
    tolTwo = -1; 
    
    maxIT1 = iterations; 
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

            % forward step
            update = (FORWARD)\(F*I+RHS); 
            update = update.*interactionMap;
            update = (1-w)*I + w*update;
            
            
            I = update;
            
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

        
       I = bsxfun(@rdivide, I, sum(I, 2));
        
        diffOne = full(max(max(abs(prevOne - I))));
        if (diffOne < tolOne)
            dispif(true, '\nrepPureGS converged after %d sequences\n', j);
            dispif(true, 'total number of iterations: %d\n', totalIter);
            break;
        end
        
        if (j==maxIT1)
            dispif(true, '\nrepPureGS did not converge after %d sequences\n', j);
            dispif(true, 'total number of iterations: %d\n', totalIter);
        end
        
        %plotToolbar(G, I)
        %FRAME(j) = getframe;
    end
    
    %figure
    %movie(FRAME, 1, 2)
    
    i = j; 
    I = I.*interactionMap;
    I = bsxfun(@rdivide, I, sum(I, 2));
    
    

    
        
    
    
    
    basis = struct('R', R, 'B', I, 'type', 'rsb');
end


    
