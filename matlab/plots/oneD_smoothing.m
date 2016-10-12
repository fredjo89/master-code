clc; clear all; close all; 


M = 3; 
N = 3*40; 


G = cartGrid([N 1]); 
G = computeGeometry(G);

rock.perm = logNormLayers(G.cartDims,'sigma', 5);
%rock.perm = ones(G.cells.num,1);


%rock.perm = temp; 

pv = partitionUI(G, [M, 1]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
%CG = storeInteractionRegionCart(CG); 
CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);
hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);
%A = A./max(max(A)); 


fSIZE = 15; 
fSIZE_2 = 12.5;
U = 2; 

%FigHandle = figure('Position', [200, 1000, 1000, 600]);
FigHandle = figure('Position', [200, 1000, 20*29, 10*29]);
area([1:N],rock.perm)
%LEG1 = legend('','');
%xlabel('');
%ylabel('');
axis([1,N,0,max(rock.perm)]);
xlabel('Cells','FontSize',fSIZE)
ylabel('Permeability','FontSize',fSIZE)
set(gca,'FontSize',fSIZE_2)
%set(gca,'xtick',x);
%set(gca,'ytick',[0:8:48]);
%set(LEG1);


basis1 = getJBasis(CG, A, 0,-1, 2/3);
basis2 = getJBasis(CG, A, 1,-1, 2/3);
basis3 = getJBasis(CG, A, 4,-1, 2/3);
basis4 = getJBasis(CG, A, 10,-1, 2/3);
basis5 = getJBasis(CG, A, 20,-1, 2/3);
basis6 = getJBasis(CG, A, 50,-1, 2/3);
basis7 = getJBasis(CG, A, 100,-1, 2/3);
basis8 = getJBasis(CG, A, 100000,-1, 2/3);



support = zeros(G.cells.num,1);

for i = 21:100
    support(i) = 1; 
end

color_1 = [120 122 112] ./ 255;

%FigHandle = figure('Position', [1200, 1000, 1000, 600]);
FigHandle = figure('Position', [1200, 1000, 20*29, 11.5*29]);
%plot(basis1.B(:,2), 'x--')
hold on
stairs([0:N-1], support,'--','color',color_1, 'linewidth',U)
plot([1:N],basis1(:,2),'k','linewidth',U)
%plot([1:N],basis2(:,2),'linewidth',U)
plot([1:N],basis3(:,2),'linewidth',U)
%plot([1:N],basis4(:,2),'linewidth',U)
plot([1:N],basis5(:,2),'linewidth',U)
%plot([1:N],basis6(:,2),'linewidth',U)
plot([1:N],basis7(:,2),'linewidth',U)
plot([1:N],basis8(:,2),'linewidth',U)
set(gca,'FontSize',fSIZE_2)

axis([1,120,0,max(basis1(:,2))])
LEG1 = legend('support', 'initial', '4 iterations', '20 iterations', '100 iterations', 'converged');

xlabel('Cells','FontSize',fSIZE);
ylabel('Basis value','FontSize',fSIZE);
set(LEG1,'FontSize',9);








AP_1 = abs(full(A*basis1(:,2)));
AP_2 = abs(full(A*basis2(:,2)));
AP_3 = abs(full(A*basis3(:,2)));
AP_4 = abs(full(A*basis4(:,2)));
AP_5 = abs(full(A*basis5(:,2)));
AP_6 = abs(full(A*basis6(:,2)));
AP_7 = abs(full(A*basis7(:,2)));
AP_8 = abs(full(A*basis8(:,2)));


U = 2; 

%FigHandle = figure('Position', [1200, 200, 1000, 600]);
FigHandle = figure('Position', [1200, 200, 15*29, 11.5*29]);
%plot(basis1.B(:,2), 'x--')
hold on
plot([1:N],AP_1, 'linewidth',U)
plot([1:N],AP_2, 'linewidth',U)
plot([1:N],AP_3, 'linewidth',U)
plot([1:N],AP_8, 'linewidth',U)
axis([1,120,0,max(AP_1)])
LEG1 = legend( 'initial', '2 iterations', '4 iterations', 'converged');
set(LEG1);
xlabel('Cells','FontSize',fSIZE);
ylabel('|AP|','FontSize',fSIZE);
set(gca,'FontSize',fSIZE_2)



