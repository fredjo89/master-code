clc; clear all; close all; 



x = [0 :0.01: 10]; 

a = 0.008; 

y = x - a*x.^3 ;



my_green_1 = [93 148 111] ./ 255;
my_green_2 = [87 160 37] ./ 255;

my_blue_1 = [61 97 209] ./ 255;
my_blue_2 = [89 89 224] ./ 255;

my_red_1 = [223 95 88] ./ 255;
my_red_2 = [193 8 23] ./ 255;



plot(x,y, 'color',my_red_2, 'linewidth',3)
hold on; 
plot(x,x, '--', 'color',my_blue_1, 'linewidth',3)

axis equal;
axis([0,10,0,10]);
xlabel('Prosessor units');
ylabel('Speedup');
set(gca,'fontsize',15)
set(gca,'xtick',[])
set(gca,'ytick',[])




LEG1 = legend('Realistic speddup','Perfect speedup');
set(LEG1);













