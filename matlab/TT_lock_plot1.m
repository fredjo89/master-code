clc; clear all; close all; 

fn = '/home/shomeb/f/fredjoha/Desktop/master-results/TT_results_1.txt'

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
    
   
    a1 = (X(11,2)/X(1,2))^(1/10)
    a2 = (X(11,3)/X(1,3))^(1/10)
    a3 = (X(11,4)/X(1,4))^(1/10)
    a4 = (X(11,5)/X(1,5))^(1/10)
    a5 = (X(11,6)/X(1,6))^(1/10)
    a6 = (X(11,7)/X(1,7))^(1/10)
    a7 = (X(11,8)/X(1,8))^(1/10)
    a8 = (X(11,9)/X(1,9))^(1/10)
    
    
    
    
    log(a1)/log(a3)
    log(a1)/log(a5)
    log(a1)/log(a7)
    
    log(a2)/log(a4)
    log(a2)/log(a6)
    log(a2)/log(a8)
    
    
    
    IT1 = 26;
    IT2 = 83;
    IT3 = 156;
    L = [
        500 500 500 500 500 500;
        500 500 500 500 500 500;
        500 500 500 500 500 500;
        ];
    for i = IT1:500
        if X(IT1,2)>X(i,4) && L(1,1)==500
            L(1,1) = i; 
        end
        if X(IT1,3)>X(i,5) && L(1,2)==500
            L(1,2) = i; 
        end
        if X(IT1,2)>X(i,6) && L(1,3)==500
            L(1,3) = i; 
        end
        if X(IT1,3)>X(i,7) && L(1,4)==500
            L(1,4) = i; 
        end
        if X(IT1,2)>X(i,8) && L(1,5)==500
            L(1,5) = i; 
        end
        if X(IT1,3)>X(i,9) && L(1,6)==500
            L(1,6) = i; 
        end
        
        if X(IT2,2)>X(i,4) && L(2,1)==500
            L(2,1) = i; 
        end
        if X(IT2,3)>X(i,5) && L(2,2)==500
            L(2,2) = i; 
        end
        if X(IT2,2)>X(i,6) && L(2,3)==500
            L(2,3) = i; 
        end
        if X(IT2,3)>X(i,7) && L(2,4)==500
            L(2,4) = i; 
        end
        if X(IT2,2)>X(i,8) && L(2,5)==500
            L(2,5) = i; 
        end
        if X(IT2,3)>X(i,9) && L(2,6)==500
            L(2,6) = i; 
        end
        
        if X(IT3,2)>X(i,4) && L(3,1)==500
            L(3,1) = i; 
        end
        if X(IT3,3)>X(i,5) && L(3,2)==500
            L(3,2) = i; 
        end
        if X(IT3,2)>X(i,6) && L(3,3)==500
            L(3,3) = i; 
        end
        if X(IT3,3)>X(i,7) && L(3,4)==500
            L(3,4) = i; 
        end
        if X(IT3,2)>X(i,8) && L(3,5)==500
            L(3,5) = i; 
        end
        if X(IT3,3)>X(i,9) && L(3,6)==500
            L(3,6) = i; 
        end
    end
   
    
    
   
    FigHandle = figure('Position', [1200, 200, 500, 500]);
    
    fSize = 10; 
    legSize = 10; 
    
    
    
    subplot(2,1,1);
    hold on; 
    plot(X(:,1), log(X(:,2)), 'b', 'LineWidth', 1)
    plot(X(:,1), log(X(:,4)), 'r', 'LineWidth', 1)
    plot(X(:,1), log(X(:,6)), 'g', 'LineWidth', 1)
    plot(X(:,1), log(X(:,8)), 'c', 'LineWidth', 1)
    LEG1 = legend('k=1', 'k=2', 'k=5', 'k=10');
    xlabel('Iterations', 'FontSize',fSize)
    ylabel('log(\infty-norm)', 'FontSize',fSize)
    set(LEG1,'FontSize',legSize);
    
    subplot(2,1,2);
    hold on; 
    plot(X(:,1), log(X(:,3)), 'b', 'LineWidth', 1)
    plot(X(:,1), log(X(:,5)), 'r', 'LineWidth', 1)
    plot(X(:,1), log(X(:,7)), 'g', 'LineWidth', 1)
    plot(X(:,1), log(X(:,9)), 'c', 'LineWidth', 1)
    xlabel('Iterations', 'FontSize',fSize)
    ylabel('log(2-norm)', 'FontSize',fSize)
%{
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters lock50
%}


    
    
    