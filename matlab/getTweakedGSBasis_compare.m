function [error] = getTweakedGSBasis_compare(CG, A, iterations, w, basisSol)

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
    
    error = zeros(iterations,3); 
    
    for i = 1:iterations
        update = (FORWARD)\(F*I+RHS); 
        update = update.*interactionMap;
        update = (1-w)*I + w*update;
        I = update;
        I = bsxfun(@rdivide, I, sum(I, 2));
        
        diff = abs(basisSol - I);
        error(i,1) = max(max(diff));                        % max inf error of the difference matrix
        error(i,2) = sqrt(max(sum(diff.^2)./sum(I.^2)));    % max relative two norm of the difference matrix
        error(i,3) = max(sum(diff)./sum(I));                % max relative 1-norm of the difference matrix

    end

    

    
    

    
    
end


    
