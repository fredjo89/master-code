% Only created to be compared with getBLGSBasis.m and see if they match. 

function [I] = getBLGSBasis_component(CG, A, iterations,tol, w)
    N = size(A,1);
    

    
    % Create interaction region 
    lens = cellfun(@numel, CG.cells.interaction);
    blocks = rldecode((1:CG.cells.num)', lens);
    interaction = vertcat(CG.cells.interaction{:});
    interactionMap = sparse(interaction, blocks, ones(size(interaction)), CG.parent.cells.num, CG.cells.num);
    G = CG.parent;
    I = controlVolumeRestriction(CG.partition)';
    
    [offsets, support, celltypes] = getGridData(CG);
    support = support+1; 
    
    
    boundary = zeros(N,1); 
    
    for i = 1:length(support)
        if celltypes(i)~=0
            boundary(support(i))=1; 
        end
    end
    
    
        
    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    M = D_inv*A;
    

    
    
    M = CG.cells.num; 
    
    for it = 1:iterations
        
        for j = 1:M
            for i = 1:N
                if boundary(i)==0
                    temp = 0; 
                    for k = 1:N
                        if (k~=i)
                            temp = temp + A(i,k)*I(k,j); 
                        end
                        I(i,j) = -temp; 
                    end
                end
            end
        end
        
        prev = I; 
        
        for j = 1:M
            for i = 1:N
                if boundary(i)==1
                    temp = 0; 
                    for k = 1:N
                        if (k~=i)
                            temp = temp + A(i,k)*prev(k,j); 
                        end
                        I(i,j) = -temp; 
                    end
                end
            end
        end
        
        
        
        I = I.*interactionMap;
        
        I = bsxfun(@rdivide, I, sum(I, 2));

        

    end

end


    
