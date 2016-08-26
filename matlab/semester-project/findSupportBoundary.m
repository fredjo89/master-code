function [ boundary ] = findSupportBoundary( currentBasis, nx, ny )

currentBasis

for i=1:nx*ny
    boundary(i) = 0; 
    if ~any(abs(currentBasis-i)<1e-10)
        south =0; north =0; west = 0; east = 0; 
        if i<=nx
            south = 1; 
        elseif i>nx*ny-nx
            north = 1; 
        end
        
        if mod(i,nx)==0
            east = 1; 
        elseif mod(i,nx)==1
            west = 1; 
        end
        
        
        
        
        if (~south && ~north && ~west && ~east)
          if ( any(abs(currentBasis-i+1)<1e-10) || any(abs(currentBasis-i-1)<1e-10) ||any(abs(currentBasis-i+nx)<1e-10)||any(abs(currentBasis-i-nx)<1e-10) )
              boundary(i)=1;
          end
        elseif west && south
            if ( any(abs(currentBasis-i+1)<1e-10)  ||any(abs(currentBasis-i+nx)<1e-10) )
              boundary(i)=1;
            end
        elseif south && east
            if (  any(abs(currentBasis-i+nx)<1e-10) )
              boundary(i)=10;
          end
        elseif north && east
            if ( any(abs(currentBasis-i-1)<1e-10) ||any(abs(currentBasis-i-nx)<1e-10) )
              boundary(i)=1;
          end
        elseif north && west
            if ( any(abs(currentBasis-i-nx)<1e-10) )
              boundary(i)=10;
          end
        elseif west
            if (  any(abs(currentBasis-i+1)<1e-10) )
              boundary(i)=1;
          end
        elseif east 
            if ( any(abs(currentBasis-i-1)<1e-10) )
              boundary(i)=1;
          end
        elseif south
           if ( any(abs(currentBasis-i+1)<1e-10) || any(abs(currentBasis-i-1)<1e-10) ||any(abs(currentBasis-i+nx)<1e-10)||any(abs(currentBasis-i-nx)<1e-10) )
              boundary(i)=1;
          end
        elseif north
           if ( any(abs(currentBasis-i+1)<1e-10) || any(abs(currentBasis-i-1)<1e-10) ||any(abs(currentBasis-i+nx)<1e-10)||any(abs(currentBasis-i-nx)<1e-10) )
              boundary(i)=1;
          end
        end
    end
end
            
              
           
        
        


end

