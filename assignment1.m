%function[] = assignment1()
clear
close all

m_e = 9.1093837015e-31;%rest mass in kg
m_n = 0.26*m_e;%Effective mass of electrons

%Nominial size of the region in nm
xlimit = [0 200e-9];
ylimit = [0 100e-9];

%Question 1: Electron Modelling
%Q1: Part 1
%What is the thermal velocity?
kb=1.380649e-23;%Boltzmann's constant [m^2*kg/K*s^2]
T = 300;%Temperature in Kelvin
vth = sqrt((kb*T)/m_n);

%Q1: Part 2
%What is the mean free path?
tmn = 0.2e-12;%mean time between collisions
R = 8.314462; %Gas Constant [J/Kmol]
vrms = sqrt((3*R*T)/m_n);
MFP = vrms*tmn;
%Q1: Part 3
%Simulation of the random motion of the electrons


%Amount of particles (total or plotted[plots first X particles])
particles = 100;
plottedparticles = 7;

ts=5e-15;%time step
tr=500e-15;%runtime

%Place particles in random locations
x = zeros([1 particles]);
y = zeros([1 particles]);
for n=1:particles
    x(n) = xlimit(1) + (xlimit(2)-xlimit(1)).*rand;
    y(n) = ylimit(1) + (ylimit(2)-ylimit(1)).*rand;
end
Vx = zeros([1 particles]);
Vy = zeros([1 particles]);

%Get random normalized angles
for n = 1:particles
    angX = cos((2*rand-1)*2*pi);
    angY = sin((2*rand-1)*2*pi);
    normalized = [angX angY]./norm([angX angY]);
    Vx(n)= vth*normalized(1);
    Vy(n)= vth*normalized(2);
end

%Physics of particles
n=1;
Px = zeros([100 particles]);
Py = zeros([100 particles]);
vAverage = 0;
for t = 0:ts:tr
    
    %Store Location
    for p = 1:particles
        Px(n, p)= x(p);
        Py(n, p)= y(p);
    end
    n=n+1;%New column for next timestep
    
    %Advance Location
    distance = Vx*ts;
    x= x + Vx*ts;
    y= y + Vy*ts;
    
    %Check boundaries
    for p = 1:particles
        if y(p) < ylimit(1) || y(p) > ylimit(2)
            Vy(p) = -1*Vy(p);
        end
        if x(p) < xlimit(1) 
            x(p) = x(p) + xlimit(2);
        elseif x(p) > xlimit(2)
            x(p) = x(p) - xlimit(2);
        end
    end
    Temp = (2/3)*(1/kb)*(1/2)*m_n*((sum(Vx)/p)^2+(sum(Vy)/p)^2);
    plot(t,Temp);
    %Temperature
    %vAverage = (vAverage + sqrt(Vx^2+Vy^2))/p;
    %Temp=(2/3)*(1/kb)*(1/2)*m_n*vAverage^2;
end

s = size(Px);
for N = 1:s(1)
    %Plot Trajectories
    for pl = 1:plottedparticles
        plot(Px(1:N, pl),Py(1:N, pl),'.');
        hold on
    end
    hold off
    %Temperature
    %vAverage = (vAverage + sqrt(Vx^2+Vy^2))/p;
    %Temp=(2/3)*(1/kb)*(1/2)*m_n*vAverage^2;
    xlim(xlimit)
    ylim(ylimit)
    xlabel('x');
    ylabel('y');
    grid on
    pause(0.01)
end
%end