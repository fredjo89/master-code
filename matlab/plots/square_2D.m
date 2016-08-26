clc; clear all; close all; 

nx = 5; 
ny = 5; 
nz = 5; 

Nx = 2; 
Ny = 2; 
Nz = 2; 

G = cartGrid([nx, ny, nz]);
rock = makeRock(G, 1, 1);
p = partitionUI(G, [Nx, Ny, Nz]);
G = mcomputeGeometry(G);


T = computeTrans(G, rock);
A = getIncomp1PhMatrix(G, T);
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);
CG = setupMexInteractionMapping(CG);

[offsets, support, types, I_com, fine, coarse] = getGridData(CG);
support = support+1; 
offsets = offsets+1; 


boundary = zeros(nx*ny*nz,1);

for i = 1:CG.cells.num
    for j = offsets(i):offsets(i+1)-1
        if types(j)==2
            boundary(support(j))=1;
        end
            
    end
end

n_tt = 0; 

for i = 1:nx*ny*nz
    if boundary(i)==1
        n_tt=n_tt+1;
    end
end



n_tt/G.cells.num

plotCellData(G,boundary); 
view(3);


figure();
slice( reshape(boundary,nx,ny,nz), [1 ny], nx, [1 nz]);
shading flat, axis equal off, set(gca, 'zdir' ,  'reverse' ), box on;
colorbar( 'horiz' );


%{
figure;
% load SPE 10 data set
mrstModule add spe10;
rock = SPE10_rock(); p=rock.poro; K=rock.perm;
% show p
slice( reshape(p,60,220,85), [1 220], 60, [1 85]);
shading flat, axis equal off, set(gca, ' zdir' ,  'reverse' ), box on;
colorbar( 'horiz' );
%}











