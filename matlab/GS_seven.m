%% 1D example 
clc; clear all; close all; 


mrstModule add spe10
layers = 35:35;
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





iterations = 1000; 
tol = 10^(-6); 
w = 0.95;

Jbasis_CONV = getJBasis(CG, A, 1000000, tol, w);

Jerror = getJBasis_compare(CG, A, iterations,  w,Jbasis_CONV.B );
RBerror = GSbasisRedBlack2(CG, A, iterations, -1, 1, Jbasis_CONV.B);


temp2 = zeros(4,iterations);
for i = 2:iterations+1
    temp2(1,i-1) = (Jerror(1,i)/Jerror(1,1))^(1/(i-1));
    temp2(2,i-1) = (Jerror(2,i)/Jerror(2,1))^(1/(i-1));
    temp2(3,i-1) = (RBerror(1,i)/RBerror(1,1))^(1/(i-1));
    temp2(4,i-1) = (RBerror(2,i)/RBerror(2,1))^(1/(i-1));
end

format long; 
temp2(1,iterations-1)
%temp2(2,iterations-1)
temp2(3,iterations-1)
%temp2(4,iterations-1)

log(temp2(1,iterations))/log(temp2(3,iterations))


x = [0:iterations]; 

lineSize = 2; 
legSize = 15; 
axisSize = 10; 
figure();
semilogy(x,Jerror(1,:),'b', 'LineWidth',lineSize);
hold on
semilogy(x,RBerror(1,:),'r', 'LineWidth', lineSize);
LEG1 = legend('Jacobi', ' Red-black');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);

%{
figure();
semilogy(x,Jerror(2,:),'b', 'LineWidth',lineSize);
hold on
semilogy(x,RBerror(2,:),'r', 'LineWidth', lineSize);
LEG1 = legend('Jacobi', ' Red-black');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Error', 'FontSize',legSize);
%}


x = [ 1:iterations]; 
figure();
semilogx(x,temp2(1,:),'b', 'LineWidth',lineSize);
hold on
semilogx(x,temp2(3,:),'r', 'LineWidth', lineSize);
LEG1 = legend('Jacobi', ' Red-black');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Convergence factor', 'FontSize',legSize);

%{
figure();
plot(temp2(2,:),'b', 'LineWidth',lineSize);
hold on
plot(temp2(4,:),'r', 'LineWidth', lineSize);
LEG1 = legend('Jacobi', ' Red-black');
set(gca, 'FontSize',axisSize);
set(LEG1,'FontSize',legSize);
xlabel('Iterations', 'FontSize',legSize);
ylabel('Convergence factor', 'FontSize',legSize);
%}



Jerror;
RBerror;

x = [0:iterations];
yolo = 0;
for i = 2:iterations+1
    temp = -1;
    for j = 1:iterations+1
       if RBerror(2,j)<=Jerror(2,i) 
           temp = j;
           break;
       end
    end
    yolo(i) = j-1;
end

plot(x,yolo,'*')
hold on
plot(x,x)
for i = 2:iterations+1
    yolo(i) = yolo(i)/x(i); 
end

yolo(1) = 1; 
figure(); 
semilogx(x,yolo)





close all;
plotCellData(G,full(Jbasis_CONV.B(:,8)))
shading flat, axis equal off, set(gca, 'zdir' ,' reverse'), box on;
colorbar( 'horiz' );
outlineCoarseGrid(G,pv)

%{
% load SPE 10 data set
mrstModule add spe10;
rock = SPE10_rock(); p=rock.poro; K=rock.perm;
% show p
slice( reshape(p,60,220,85), [1 220], 60, [1 85]);
shading flat, axis equal off, set(gca, 'zdir' ,' reverse'), box on;
colorbar( 'horiz' );
% show Kx
slice( reshape(log10(K(:,1)),60,220,85), [1 220], 60, [1 85]);
shading flat, axis equal off, set(gca, 'zdir' , 'reverse' ), box on;
h=colorbar( 'horiz' );
set(h, 'XTickLabel' ,10.^[get(h, 'XTick' )]);
set(h, 'YTick' ,mean(get(h, 'YLim' )), 'YTickLabel' , 'mD' );
%}




