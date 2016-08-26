function [I, error] = getJBasis(CG, A, iterations,tol, w)

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
   
    error_init = full(sum(sum(abs(A*I))));
    error(1) =  1;
    
    tic
    
    for  i =1:iterations
        
        prev = I; 
        
        update = -M*I;
        update = update.*interactionMap;
        
        I = I + w*update;
        I = bsxfun(@rdivide, I, sum(I, 2));
        
        error(i+1) = sum(sum(abs(A*I)))/error_init;
        
        if abs(error(i+1)-error(i))<tol
            X = ['TOL, J: ', num2str(i)];
            disp(X);
            break;
        elseif error(i+1)>error(i)
            X = ['INCREASE, J: ', num2str(i)];
            disp(X);
            break;
        end
        
        if rem(i,1000)==0
            disp([num2str(i),', errordiff: ',num2str(abs(error(i+1)-error(i))) ]);
            toc
        end
    end


    
end


    
