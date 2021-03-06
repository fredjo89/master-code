function [I, error] = getBLGSBasis(CG, A, iterations,tol, w)
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
   
    
   RHS = repmat(boundary, 1,size(I,2));
    
   
   %{
   disp('yolo')
   tic
   temp = repmat(boundary, 1,N);
   toc
   tic
   A2 = A.*(1-temp);
   toc
   tic
   temp2 = spdiags(boundary, 0, N, N);
   toc
   tic
   A2=A2+temp2;
   toc
   disp('yolo')
   
   clear temp; 
   clear temp2; 
   %}
   
   boundary = full(boundary);
   
   temp = zeros(1,N);
   tic
    A2 = A; 
    for i = 1:N
        if boundary(i)==1
            A2(i,:) = temp;
            A2(i,i) = 1; 
        end
        if rem(i,1000)==0
            i
            toc
        end
        
    end
  toc
  
  clear temp;
    
   

    D = diag(A);
    n = numel(D);
    D_inv = spdiags(1./D, 0, n, n);
    M = D_inv*A;
    
    D2 = diag(A2);
    F2 = - (triu(A2) - spdiags(D2, 0, n, n));
    FORWARD = A2+F2;
    
    error_init = full(sum(sum(abs(A*I))));
    error(1) =  1;
    
    tic
    for i = 1:iterations
        prev = I; 
        
        temp1 = (FORWARD)\(F2*I+RHS.*I);
        
        if (full(sum(sum((temp1-I).*RHS)))~=0)
            disp(['Something went wrong']);
            break; 
        end
        
        temp3 = -M*temp1; 
        temp3 = temp3.*RHS; 
        
        I = temp1 + temp3; 
        
        I = I.*interactionMap;
        I = bsxfun(@rdivide, I, sum(I, 2));
        
      
        error(i+1) = sum(sum(abs(A*I)))/error_init;
        
        if abs(error(i+1)-error(i))<tol
            X = ['TOL, BLGS: ', num2str(i)];
            disp(X);
            break;
        elseif error(i+1)>error(i)
            X = ['INCREASE, BLGS: ', num2str(i)];
            disp(X);
            break;
        end
        
        if rem(i,1000)==0
            toc
            disp([num2str(i),', errordiff: ',num2str(abs(error(i+1)-error(i))) ]);
        end
        
    end


end


    
