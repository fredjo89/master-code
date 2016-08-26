clc; clear all; close all; 


N = 25;                            % Number of cells
n = 5;                             % Number of cells in each block
M = ceil(N / n);                    % Number of blocks
G = cartGrid([N N]); 
G = computeGeometry(G);
rock = makeRock(G, 1, 1);
pv = partitionUI(G, [M, M]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);
hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);



%{
mrstModule add spe10
layers = 85:85;
[G, ~, rock] = SPE10_setup(layers);
pv = partitionUI(G, [6, 11, ceil(G.cartDims(3)./5)]);
CG = generateCoarseGrid(G,pv); 
CG = coarsenGeometry(CG); 
CG = storeInteractionRegionCart(CG);
%CG = storeInteractionRegionCart(CG, 'adjustCenters', false, 'edgeBoundaryCenters', false);
CG = setupMexInteractionMapping(CG);
hT = computeTrans(G,rock); 
A = getIncomp1PhMatrix(G, hT);

plotCellData(G,log10(rock.perm(:,1)))
shading flat, axis equal off, set(gca, 'zdir' ,' reverse'), box on;
colorbar( 'horiz' );
outlineCoarseGrid(G,pv)
%}




iterations = 1000000; 

Jerror1 = getJBasis_compare2(CG, A, iterations,  2/3);
Jerror2 = getJBasis_compare2(CG, A, iterations,  0.95);
RBerror = RedBlack_compare(CG, A, iterations, 1);

AST = min([ min(Jerror1) min(Jerror2) min(RBerror) ]);

disp(['Start Error: ', num2str(Jerror1(1))]);
disp(['AST: ', num2str(AST)]);
disp(['Jerror1: ', num2str(min(Jerror1)-AST)]);
disp(['Jerror2: ', num2str(min(Jerror2)-AST)]);
disp(['RBerror: ', num2str(min(RBerror)-AST)]);

Jerror1 = (Jerror1-AST)/(Jerror1(1)-AST);
Jerror2 = (Jerror2-AST)/(Jerror2(1)-AST);
RBerror = (RBerror-AST)/(RBerror(1)-AST);


F = 1;
for i = 2:length(Jerror2)
    temp = -1;
    for j = 1:length(RBerror)
       if RBerror(j)<=Jerror2(i) 
           temp = j;
           break;
       end
    end
    F(i) = j-1;
end

for i = 2:length(Jerror2)
    F(i) = F(i)/(i-1); 
end

L1 = length(Jerror1);
L2 = length(Jerror2);
L3 = length(RBerror);





lineSize = 2; 
legSize = 15; 
axisSize = 10; 
figure();
R3 = 50000; 
R1 = R3 + L1-L3; 
R2 = R3 + L2 - L3;
semilogy([0:L3-1 - R3],RBerror(1:L3-R3),'k','LineWidth',lineSize)
hold on
semilogy([0:L1-1-R1],Jerror1(1:L1-R1),'r','LineWidth',lineSize)
semilogy([0:L2-1-R2],Jerror2(1:L2-R2),'b','LineWidth',lineSize)
LEG1 = legend('Red-black','Jacobi, 2/3','Jacobi, 0.95');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);





L = length(RBerror)-1; 
start = 3; 
lineSize = 2; 
legSize = 15; 
axisSize = 10; 
figure();
REMOVE = 200;
semilogx([start:L-REMOVE-1],F(start+1:L-REMOVE),'-','LineWidth',lineSize)
set(gca, 'FontSize',axisSize);
%set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Fraction', 'FontSize',legSize);








%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters error_spe85
%}

%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters fraction_400
%}


