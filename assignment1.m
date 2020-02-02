%% Beginning of Assignment 1
clear
close all

part3 = true;

m_e = 9.1093837015e-31;%rest mass in kg
m_n = 0.26*m_e;%Effective mass of electrons

%Nominial size of the region in nm
xlimit = [0 200e-9];
ylimit = [0 100e-9];
%% Part 1: Electron Modelling
%% P1: Question 1
%What is the thermal velocity?
kb=1.380649e-23;%Boltzmann's constant [m^2*kg/K*s^2]
T = 300;%Temperature in Kelvin
vth = sqrt((3*kb*T)/m_n);

%% P1: Question 2
%What is the mean free path?
tmn = 0.2e-12;%mean time between collisions
MFP = vth*tmn;


%% P1: Question 3
%Setting initial values of particles


%Amount of particles (total or plotted[plots first X particles])
particles = 1000;
plottedparticles = 7;

timesteps = 100;
ts=2e-15;%time step
tr=timesteps*ts;%runtime

%Place particles in random locations
x = zeros([1 particles]);
y = zeros([1 particles]);
for n=1:particles
    x(n) = xlimit(1) + (xlimit(2)-xlimit(1)).*rand;
    y(n) = ylimit(1) + (ylimit(2)-ylimit(1)).*rand;
    if part3
        while (x(n) >= 0.8e-7 && x(n) <= 1.2e-7 && y(n) <= 0.4e-7)...
                || (x(n) >= 0.8e-7 && x(n) <= 1.2e-7 && y(n) >= 0.6e-7)
            x(n) = xlimit(1) + (xlimit(2)-xlimit(1)).*rand;
            y(n) = ylimit(1) + (ylimit(2)-ylimit(1)).*rand;
        end
    end
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
%% P2: Question 1: Get randomized velocities from maxwell-boltzmann distribution
Fx = sqrt((m_n)/(2*pi*kb*T))*exp(((-m_n*Vx.^2)/(2*kb*T)));
Fy = sqrt((m_n)/(2*pi*kb*T))*exp(((-m_n*Vy.^2)/(2*kb*T)));
% Vx = Fx;
% Vy = Fy;
%% P1: Question 3
%%Simulation of the random motion of the electrons after time t

%Initialize variables used for storage for plotting later
n=1;
Px = zeros([timesteps particles]);
Py = zeros([timesteps particles]);
Temp = zeros([1 timesteps]);
Time = zeros([1 timesteps]);
for t = 0:ts:tr
    
    %Store Location
    for p = 1:particles
        Px(n, p)= x(p);
        Py(n, p)= y(p);
    end
    
    %Advance Location
    distance = Vx*ts;
    x= x + Vx*ts;
    y= y + Vy*ts;
    
    %Check boundaries
    for p = 1:particles
        if y(p) <= ylimit(1) || y(p) >= ylimit(2)
            Vy(p) = -1*Vy(p);
        end
        if x(p) <= xlimit(1) 
            x(p) = x(p) + xlimit(2);
        elseif x(p) >= xlimit(2)
            x(p) = x(p) - xlimit(2);
        end
    end
    
    %Check Box Boundaries
    if part3
         %% Part 3: Enhancements Question 1 Boundary boxes
         for p = 1:particles
            if (x(p) >= 0.8e-7 && y(p) <= 0.4e-7 && x(p) <= 0.9e-7)... 
                && (x(p) >= 0.8e-7 && y(p) >= 0.6e-7 && x(p) <= 0.9e-7)...
                && (x(p) <= 1.2e-7 && y(p) <= 0.4e-7 && x(p) <= 1.1e-7)...
                && (x(p) <= 1.2e-7 && y(p) >= 0.6e-7 && x(p) <= 1.1e-7)
                Vx(p) = -1*Vx(p);
            elseif (x(p) >= 0.8e-7 && x(p) <= 1.2e-7 && y(p) <= 0.4e-7 && y(p) >= 0.35e-7)...
                    && (x(p) >= 0.8e-7 && x(p) <= 1.2e-7 && y(p) >= 0.6e-7 && y(p) <= 0.65e-7)
                Vy(p) = -1*Vy(p);
            end
        end
    end
        
        
    Temp(n) = (2/3)*(1/kb)*(1/2)*m_n*(sum(Vx.^2+Vy.^2))/p;
    Time(n) = t;
    n=n+1;%New column in matrix for each timestep
end

%% P1: Question 3 i) Plot Trajectories
s = size(Px);
figure;
for N = 1:s(1)
    for pl = 1:plottedparticles
        plot(Px(1:N, pl),Py(1:N, pl),'.');
        hold on
    end
    hold off
    xlim(xlimit)
    ylim(ylimit)
    xlabel('x');
    ylabel('y');
    grid on
    if part3
        %% Part 3: Enhancements Question 1 Boundary boxes
        rectangle('Position',[0.8e-7 0 0.4e-7 0.4e-7])
        rectangle('Position',[0.8e-7 0.6e-7 0.4e-7 0.4e-7])
    end
    pause(0.001)
end

%% P1: Question 3 ii) Temperature plot
figure;
plot(Time,Temp,'-')
xlabel('Time (s)');
ylabel('Temperature (K)');
title('Temperature Plot')
grid on

%% Part 2: Collisions with Mean Free Path
%Velocity from part 1 is now assigned random velocities using
%Maxwell-Boltzmann, average speed should be vth
figure;
histogram(Fx,20)
figure;
histogram(Fy,20)


%% Part 3: Enhancements


%% End of assignment