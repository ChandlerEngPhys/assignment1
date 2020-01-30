%function[] = assignment1()
clear
clf
%Effective mass of electrons
m_n = 0.26;

%Nominial size of the region in nm
xlimit = [0 200];
ylimit = [0 100];

%question1();

%end

%function[] = question1()
%Question 1 Electron Modelling
%Part 1
%What is the thermal velocity?
%Temperature in Kelvin
T = 300;
MFP = 1;

%Part 2
%What is the mean time between collisions?

%Part 3
grid ON

particles = 5;
V0= 1;

ts=1;%time step
tr=100;%runtime

% xlim(axes, xlimit);
% ylim(axes, ylimit);
%plot(x,y,'o');
%Vx(1:nTraj) = V0 * cos(InitalAngle);
%Vy(1:nTraj) = V0 * sin(InitalAngle);
%Vx = Vx + dvx;
%dx = Vx * dt;
x = zeros([1 particles]);
y = zeros([1 particles]);
for n=1:particles
    x(n) = xlimit(1) + (xlimit(2)-xlimit(1)).*rand;
    y(n) = ylimit(1) + (ylimit(2)-ylimit(1)).*rand;
end
Vx = zeros([1 particles]);
Vy = zeros([1 particles]);

for n = 1:particles
    Vx(n)= V0*cos((2*rand-1)*2*pi);
    Vy(n)= V0*sin((2*rand-1)*2*pi);
end

n=1;
Px = zeros([tr/ts particles]);
Py = zeros([tr/ts particles]);
for t = 1:ts:tr
    
    %Store Location
    for p = 1:particles
        Px(n, p)= x(p);
        Py(n, p)= y(p);
    end
    
    %Advance Location
    x= x + Vx*ts;
    y= y + Vy*ts;
    
    %Check boundaries
    for p = 1:particles
        if y(p) < ylimit(1) || y(p) > ylimit(2)
            Vy(p) = -1*Vy(p);
        end
        if x(p) < xlimit(1) 
            x(p) = xlimit(2);
        elseif x(p) > xlimit(2)
            x(p) = xlimit(1);
        end
    end
    n=n+1;
end
for N = 1:length(Px)
    %Plot Location
    for p = 1:particles
%         hold on
%         plot(Px(1:N, p),Py(1:N, p));
%         hold off
        plot(Px(1:N, 1),Py(1:N, 1));
        hold on
        plot(Px(1:N, 2),Py(1:N, 2));
        hold off
    end
%     plot(Px(1:N, 1),Py(1:N, 1));
%     hold on
%     plot(Px(1:N, 2),Py(1:N, 2));
%     hold off
    
    xlabel('x');
    ylabel('y');
    xlim(xlimit)
    ylim(ylimit)
    grid on
    pause(0.01)
end
%end