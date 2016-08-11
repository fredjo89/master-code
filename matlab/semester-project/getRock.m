function [ rock ] = getRock( G, nCellsX, nCellsY,midWidth, channelWidth,...
                    outerPoro, innerPoro, channelPoro )

% Create Permeability and porosity fields
p = gaussianField(G.cartDims, outerPoro, [3 3 3] , 2.5);
K = p.^3*10^(-10)./(0.81*72*(1-p).^2);
rock.poro = p(:);
rock.perm = K(:);
nCellsXHalf = ceil(nCellsX/2);
jStart = 1 - ceil(midWidth/2); 
jEnd   = jStart + midWidth-1; 
for j = jStart:jEnd
    p = gaussianField([ 1 nCellsY], innerPoro, [3 3 3] , 2.5);
    p = p(:);
    K = p.^3*10^(-10)./(0.81*72*(1-p).^2);
    for i = 1:nCellsX
        rock.poro(i*nCellsY -nCellsXHalf+j+1) = p(i);  
        rock.perm(i*nCellsY -nCellsXHalf+j+1) = K(i);
    end
end
kStart = 1 - ceil(channelWidth/2)-1; 
kEnd   = kStart + channelWidth-1;
for k = kStart:kEnd
    p = gaussianField([ 1 midWidth], channelPoro, [3 3 3] , 2.5);
    p = p(:);
    K = p.^3*10^(-10)./(0.81*72*(1-p).^2);
    
    for j = jStart:jEnd
        rock.poro((nCellsXHalf+k)*nCellsX+nCellsXHalf+j) = p(j+1-jStart);
        rock.perm((nCellsXHalf+k)*nCellsX+nCellsXHalf+j) = K(j+1-jStart);
    end
    
end


end

