function [I, error] = getRedBlackBasis(CG, A, iterations,tol, w)

    % Create interaction region 
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);
    G = CG.parent;
    I = controlVolumeRestriction(CG.partition)';
    
    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    M = D_inv*A;
    
    [offsets, support, celltypes] = getGridData(CG);

    support = support+1; 
    
    N = G.cells.num;
    temp = M; 
    temp(temp~=0) = 1;
    temp = temp - speye(N);
    [col,row] = find(temp');
    
    B = zeros(1,N);
    for i = 1:length(support)
        if (celltypes(i)==2)
            B(support(i))=1; 
        end
    end
    
    l = length(col); 
    k = 1; 
    for i = 1 : N
        if (B(i)==0)
            if (row(k)>i)
                break; 
            end
            while (row(k)<i && k<l)
                k=k+1;
            end
            while (row(k) ==i && k<l)
                B(col(k))=1;
                k = k+1;
            end
        end 
    end
    B = B';
    BLACK = B; 
    for i = 2:CG.cells.num
        BLACK = [BLACK B];
    end
    RED = ones(N,CG.cells.num) - BLACK;
    
    error_init = full(sum(sum(abs(A*I))));
    error(1) =  1;
    
    for i = 1:iterations
        
        prev = I; 
        
        % Red cells 
        update = -M*I; 
        update = update.*RED;
        I = I + w*update;
        
        % Black cells 
        update = -M*I; 
        update = update.*interactionMap;
        update = update.*BLACK; 
        I = I + w*update;
        I = bsxfun(@rdivide, I, sum(I, 2));
        
        error(i+1) = sum(sum(abs(A*I)))/error_init;
        
        if abs(error(i+1)-error(i))<tol
            X = ['TOL, RB: ', num2str(i)];
            disp(X);
            break;
        elseif error(i+1)>error(i)
            X = ['INCREASE, RB: ', num2str(i)];
            disp(X);
            break;
        end
        
        if rem(i,1000)==0
            disp([num2str(i),', errordiff: ',num2str(abs(error(i+1)-error(i))) ]);
        end
        
    end

    

end


    
