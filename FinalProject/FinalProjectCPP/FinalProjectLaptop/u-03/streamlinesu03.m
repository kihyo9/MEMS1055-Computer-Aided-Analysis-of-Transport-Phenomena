data = readtable('info-run-28-04-2020 03-30-35.txt');
matrix = data{:,:};
%matrix(94,4) = matrix(94,4) - 4;
count = linspace(1,1157,1157);

semilogy(count, matrix(:,1), count, matrix(:,2), count, matrix(:,3), count, matrix(:,4))
legend('u','v','stream','vorticity','Location','southwest');
xlabel('Timestep #')
ylabel('Average change per gridpoint')

u = readtable('u-fix.txt');
u = u{:,:};
u = u';
u = flip(u,1);
v = readtable('v-fix.txt');
v = v{:,:};
v = v';
v = flip(v,1);
x = linspace(0,9,91);
y = linspace(0,1,11);

starty = linspace(0,1,31);
startx = ones(31,1)*1.6;

figure()
hold on
quiver(x,y,u,v)
h=fill([0,1.5,1.5,0],[0,0,0.5,0.5],'black');
hold off
figure()
h=fill([0,1.5,1.5,0],[0,0,0.5,0.5],'black');
streamline(x,y,u,v,startx,starty)
% starty = linspace(0,0.3,4);
% startx = ones(4,1)*1.5;
% streamline(x,y,u,v,startx,starty)
