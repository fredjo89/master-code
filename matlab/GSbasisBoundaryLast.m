function [error, i] = GSbasisBoundaryLast(CG, A, iterations, w)
    A = A - diag(sum(A, 2));

    % Create interaction region 
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);
    G = CG.parent;
    I = controlVolumeRestriction(CG.partition)';
    
    
    [offsets, support, celltypes] = getGridData(CG);

    support = support+1; 
    
    N = size(A,1);
    
    boundary = zeros(N,1); 
    
    
    for i = 1:length(support)
        if celltypes(i)~=0
            boundary(support(i))=1; 
        end
    end
    
    RHS = boundary; 
    for i  = 2 : size(I,2)
        RHS = [RHS boundary];
    end
    
    
    
    
    A1 = A; 
    A2 = A; 
    
    
    for i = 1:N
        if boundary(i)==1
            for (j = 1:N)
                A1(i,j) = 0;
            end
            A1(i,i) = 1; 
        else
            for (j = 1:N)
                A2(i,j) = 0;
            end
            A2(i,i) = 1; 
        end
            
    end
    
    D1 = diag(A1);
    n = numel(D1);
    D1_inv = spdiags(1./D1, 0, n, n);
    F1 = - (triu(A1) - spdiags(D1, 0, n, n));
    E1 = - (tril(A1) - spdiags(D1, 0, n, n)) ;
    D1 = spdiags(D1, 0, n, n);
    FORWARD1 = A1+F1;
    BACKWARD1 = A1+E1;
    
    D2 = diag(A2);
    n = numel(D2);
    D2_inv = spdiags(1./D2, 0, n, n);
    F2 = - (triu(A2) - spdiags(D2, 0, n, n));
    E2 = - (tril(A2) - spdiags(D2, 0, n, n)) ;
    D2 = spdiags(D2, 0, n, n);
    FORWARD2 = A2+F2;
    BACKWARD2 = A2+E2;
    
   
    
    error_init = full(sum(sum(abs(A*I))));
    error(1) =  1;
    
    full(RHS)
    
    
    for i = 1:iterations
        
        
        
        I = (FORWARD1)\(F1*I+RHS.*I);
        %disp(['Distance from normalization: ', num2str(sum(sum(I')-1))]);
        
        
        I = (FORWARD2)\(F2*I+(1-RHS).*I);
   
        
        I = I.*interactionMap;
        I = bsxfun(@rdivide, I, sum(I, 2));
        
        
        
        error(i+1) = sum(sum(abs(A*I)))/error_init;
        
        if abs(error(i+1)-error(i))<10^(-12)
            X = ['TOL, GS: ', num2str(i)];
            disp(X);
            break;
        elseif error(i+1)>error(i)
            X = ['INCREASE, GS: ', num2str(i)];
            disp(X);
            break;
        end
        
        if rem(i,1000)==0
            disp([num2str(i),', errordiff: ',num2str(abs(error(i+1)-error(i))) ]);
        end
    end

end


    
