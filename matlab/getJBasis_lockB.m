function [I, error] = getJBasis_lockB(CG, A, iterations,tol, w, stop)

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
   
    % Finding boundary cells and stores them in RHS
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
    
    
    
    error_init = full(sum(sum(abs(A*I))));
    error(1) =  1;
    
    tic
    for  i =1:iterations
        prev = I; 
        
        for s = 1 : 1
            update = -M*I;
            update = update.*interactionMap; 
            update = update.*(1-RHS); 
            I = I + w*update; 
        end
        
        
        
        update = -M*I;
        update = update.*interactionMap;
        
        I = I + w*update;
        I = bsxfun(@rdivide, I, sum(I, 2));
        
        error(i+1) = sum(sum(abs(A*I)))/error_init;
        
        if (error(i+1) < stop)
            X = ['reached stop in: ', num2str(i)];
            disp(X);
            break;
        end
        
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


    
