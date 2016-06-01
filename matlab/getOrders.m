function [Orders] = getOrders(CG, A)
    
    [rows, cols ] = find (A); 
    G = CG.parent;
    
    % find dependence pattern
    currentCol = 1; 
    index = 0; 
    for (i = 1:length(rows))
        if ( cols(i)~=currentCol )
            currentCol= currentCol+1;
            index = 0; 
        end
        if (rows(i)~=currentCol)
            index = index + 1; 
            J_index(currentCol, index) = rows(i); 
        end
    end


    %find center cell
    for (it = 1:CG.cells.num)
        cCell = CG.cells.centers(it);
        order = cCell; 

        currentLength = 1;
        prevLength = 0; 
        placed = 0; 
        for (i = 1:100)
            startPlace = placed + 1;  
            endPlace = currentLength; 
            temp = [J_index(order(startPlace), :)];
            for (j = startPlace+1:endPlace)
                temp = [temp J_index(order(j), :)];
            end
            temp(temp==0) = [];
            order = [order temp];

            [Xs, SortVec] = sort(order(:));
            UV(SortVec) = ([1; diff(Xs)] ~= 0);
            UV = UV(1:length(SortVec));
            order = order(UV);

            placed = placed + 1+endPlace - startPlace; 
            currentLength = length(order);
            if currentLength == G.cells.num
                break;
            end
        end
        if (it==1)
            Orders = order'; 
        else
            Orders = [Orders order']; 
        end
    end
    


end



