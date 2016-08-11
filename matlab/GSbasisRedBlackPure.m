function [basis, i] = GSbasisRedBlackPure(CG, A, iterations, tol, w)
    A = A - diag(sum(A, 2));
    R = controlVolumeRestriction(CG.partition);
    % Create interaction region 
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);
    G = CG.parent;
    I = controlVolumeRestriction(CG.partition)';

    RHS = sparse(G.cells.num, CG.cells.num);
    for k = 1:CG.cells.num
        temp = zeros(1,G.cells.num); 
        temp(CG.cells.centers(k)) = 1; 
        A(CG.cells.centers(k),:) =temp; 
        RHS(CG.cells.centers(k),k)=1;
    end
    
    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    F = - (triu(A) - spdiags(D, 0, n, n));
    FORWARD = A+F;
    M = D_inv*A;
    
    [offsets, support, celltypes] = getGridData(CG);

    
    support = support +1; 
    
    N = G.cells.num;
    temp = M; 
    temp(temp~=0) = 1;
    temp = temp - speye(N);
    [col,row] = find(temp');
    
    B = zeros(1,N);
    
    n = sqrt(N);
    
    if (rem(n,2)==1)
       for i=2:2:N
           B(i) = 1;
       end
    else
        temp = -1; 
        for i=1:N
            if rem(i-1,n)==0
                temp = temp *-1;
            end
            
            if (temp==-1 && rem(i,2)==0)
                B(i)=1; 
            elseif(temp==1 && rem(i,2)==1)
                B(i) = 1; 
            end
        end
    end
    
    

    

  
    B = B';
    BLACK = B; 
    for i = 2:CG.cells.num
        BLACK = [BLACK B];
    end
    RED = ones(N,CG.cells.num) - BLACK;
    
    figure()
    pv = partitionUI(G, [sqrt(CG.cells.num), sqrt(CG.cells.num)]);
    colorMatrix = [ 0 0 0; 1 0 0];
    plotCellData(G,BLACK(:,1),'EdgeColor', 'y')
    plotGrid(G,'FaceColor', 'none')
    outlineCoarseGrid(G,pv,[0.7 0.7 0.7])
    axis equal tight off; 
    colormap(colorMatrix)
    %colormap(1,1) = 1
 
    
    for i = 1:iterations
        prev = I; 
        
        % Black cells 
        update = -M*I+D_inv*RHS; 
        update = update.*interactionMap;
        update = update.*BLACK; 
        I = I + update;
        
        % Red cells 
        update = -M*I+D_inv*RHS; 
        update = update.*interactionMap;
        update = update.*RED;
        I = I + update;
        I = bsxfun(@rdivide, I, sum(I, 2));
  
        % update and normalize
        I = (1-w)*prev + w*I; 
    end
    
    basis = struct('R', R, 'B', I, 'type', 'rsb');
    
end


    
