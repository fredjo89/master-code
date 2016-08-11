clc; clear all; close all; 

fn = '/home/shomeb/f/fredjoha/Desktop/master-results/TT_results_4.txt'





    fileID = fopen(fn);
    C = textscan(fileID,'%f');
    data = C{1};
    
    
    k = 1;
    for i = 1:length(data)
        if k ==10
            k = 1;
        end
        l = ceil(i/9);
        X(l,k) = data(i);
        k = k+1; 
    end
    LENGTH = 401;
    

    format long; 
    temp = LENGTH - 400; 
    temp2 = X(LENGTH,1) - X(temp,1);
    
    
    a1 = (X(LENGTH,2)/X(temp,2))^(1/temp2)
    a2 = (X(LENGTH,3)/X(temp,3))^(1/temp2)
    a3 = (X(LENGTH,4)/X(temp,4))^(1/temp2)
    a4 = (X(LENGTH,5)/X(temp,5))^(1/temp2)
    
    a5 = (X(LENGTH,6)/X(temp,6))^(1/temp2)
    a6 = (X(LENGTH,7)/X(temp,7))^(1/temp2)
    a7 = (X(LENGTH,8)/X(temp,8))^(1/temp2)
    a8 = (X(LENGTH,9)/X(temp,9))^(1/temp2)
    
    
    format short; 
    log(a1)/log(a2)
    log(a1)/log(a3)
    log(a1)/log(a4)
    
    log(a5)/log(a6)
    log(a5)/log(a7)
    log(a5)/log(a8)
    
    
    k = 1; 
    for i = 3:2:400
        Y(k,1) = X(i,1);
        Y(k,2) = (X(i,2)/X(i-2,2))^(1/10);
        Y(k,3) = (X(i,3)/X(i-2,3))^(1/10);
        Y(k,4) = (X(i,4)/X(i-2,4))^(1/10);
        Y(k,5) = (X(i,5)/X(i-2,5))^(1/10);
        Y(k,6) = (X(i,6)/X(i-2,6))^(1/10);
        Y(k,7) = (X(i,7)/X(i-2,7))^(1/10);
        Y(k,8) = (X(i,8)/X(i-2,8))^(1/10);
        Y(k,9) = (X(i,9)/X(i-2,9))^(1/10);
        k = k+1; 
    end
    
    Y(:,3) = log(Y(:,2))./log(Y(:,3))
    Y(:,6) = log(Y(:,5))./log(Y(:,6))
    format long; 
    
    hold on; 
    %plot(Y(:,1),Y(:,2))
    plot(Y(:,1),Y(:,6))
    %plot(Y(:,1),Y(:,4))
    %plot(Y(:,1),Y(:,5))
    
    
    
    
    
    
    
    
    
    
    
    
    
    IT1 = 33;
    IT2 = 201;
    IT3 = 479;
    
    for i = 1:LENGTH
        if X(i,1) == 35
            IT1 = i;
        end
        if X(i,1) == 205
            IT2 = i ;
              end
        if X(i,1) == 480
            IT3 = i ;
        end
    end
    
    L = [
        LENGTH LENGTH LENGTH LENGTH LENGTH LENGTH;
        LENGTH LENGTH LENGTH LENGTH LENGTH LENGTH;
        LENGTH LENGTH LENGTH LENGTH LENGTH LENGTH;
        ];
    for i = IT1:LENGTH
        if X(IT1,2)>X(i,3) && L(1,1)==LENGTH
            L(1,1) = X(i,1); 
        end
        if X(IT1,2)>X(i,4) && L(1,2)==LENGTH
            L(1,2) = X(i,1);
        end
        if X(IT1,2)>X(i,5) && L(1,3)==LENGTH
            L(1,3) = X(i,1);
        end
        if X(IT1,6)>X(i,7) && L(1,4)==LENGTH
            L(1,4) = X(i,1);
        end
        if X(IT1,6)>X(i,8) && L(1,5)==LENGTH
            L(1,5) = X(i,1);
        end
        if X(IT1,6)>X(i,9) && L(1,6)==LENGTH
            L(1,6) = X(i,1);
        end
        
        
        if X(IT2,2)>X(i,3) && L(2,1)==LENGTH
            L(2,1) = X(i,1);
        end
        if X(IT2,2)>X(i,4) && L(2,2)==LENGTH
            L(2,2) = X(i,1);
        end
        if X(IT2,2)>X(i,5) && L(2,3)==LENGTH
            L(2,3) = X(i,1);
        end
        if X(IT2,6)>X(i,7) && L(2,4)==LENGTH
            L(2,4) = X(i,1);
        end
        if X(IT2,6)>X(i,8) && L(2,5)==LENGTH
            L(2,5) = X(i,1);
        end
        if X(IT2,6)>X(i,9) && L(2,6)==LENGTH
            L(2,6) = X(i,1);
        end
        
        
        if X(IT3,2)>X(i,3) && L(3,1)==LENGTH
            L(3,1) = X(i,1); 
        end
        if X(IT3,2)>X(i,4) && L(3,2)==LENGTH
            L(3,2) = X(i,1);
        end
        if X(IT3,2)>X(i,5) && L(3,3)==LENGTH
            L(3,3) = X(i,1);
        end
        if X(IT3,6)>X(i,7) && L(3,4)==LENGTH
            L(3,4) = X(i,1);
        end
        if X(IT3,6)>X(i,8) && L(3,5)==LENGTH
            L(3,5) = X(i,1);
        end
        if X(IT3,6)>X(i,9) && L(3,6)==LENGTH
            L(3,6) = X(i,1);
        end
    end
   
    
    L;
   
    FigHandle = figure('Position', [1200, 200, 500, 500]);
    
    fSize = 10; 
    legSize = 10; 
    
    
    
    subplot(2,1,1);
    hold on; 
    plot(X(:,1), log(X(:,2)), 'b', 'LineWidth', 1)
    plot(X(:,1), log(X(:,3)), 'r', 'LineWidth', 1)
    plot(X(:,1), log(X(:,4)), 'g', 'LineWidth', 1)
    plot(X(:,1), log(X(:,5)), 'c', 'LineWidth', 1)
   LEG1 = legend('k=1', 'k=2', 'k=5', 'k=10');
    xlabel('Iterations', 'FontSize',fSize)
    ylabel('log(\infty-norm)', 'FontSize',fSize)
    set(LEG1,'FontSize',legSize);
    
    subplot(2,1,2);
    hold on; 
    plot(X(:,1), log(X(:,6)), 'b', 'LineWidth', 1)
    plot(X(:,1), log(X(:,7)), 'r', 'LineWidth', 1)
    plot(X(:,1), log(X(:,8)), 'g', 'LineWidth', 1)
    plot(X(:,1), log(X(:,9)), 'c', 'LineWidth', 1)
    xlabel('Iterations', 'FontSize',fSize)
    ylabel('log(2-norm)', 'FontSize',fSize)
%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters lock200
%}


    
    
    