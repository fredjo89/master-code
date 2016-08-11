function [ G, rock, p, testcase ] = getCase()
% setup case

testcase = 'custom';

switch lower(testcase)
    case 'custom'
        temp1 = 80; 
        temp2 = 2; 
        
        [nx ny] = deal(temp1); 
        theBasisDesider = ceil(temp1/temp2);
        cgxy = ceil(nx / theBasisDesider);
        G = cartGrid([nx ny]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [cgxy, cgxy]);
    case 'non_homo'
        temp1 = 120; 
        temp2 = 4; 
        [nx ny] = deal(temp1); 
        theBasisDesider = ceil(temp1/temp2);
        G = cartGrid([nx ny], [2, 2]);
        G = computeGeometry(G);
        p = gaussianField(G.cartDims, [.1 .5], [11 3] , 4.5);
        K = p.^3*10^(-10)./(0.81*72*(1-p).^2);
        rock.poro = p(:); 
        rock.perm = K(:); 
        cgxy =ceil(nx/theBasisDesider);
        p = partitionUI(G, [cgxy, cgxy]);
    case 'simple'
        G = cartGrid([15, 10]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [3, 2]);
    case 'medium'
        G = cartGrid([80, 40]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [8, 4]);
    case 'big'
        G = cartGrid([50, 50, 20]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [5, 5, 2]);
    case 'huge'
        G = cartGrid([200, 200, 100]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [10, 10, 10]);
        G = mcomputeGeometry(G);
    case 'spe10'
        mrstModule add spe10
        layers = 1:85;
        [G, ~, rock] = SPE10_setup(layers);
        p = partitionUI(G, [6, 11, ceil(G.cartDims(3)./5)]);
    otherwise
        error();
end



end

