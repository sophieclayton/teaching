% MATLAB script for OCN430 practical, October 2014
% Sophie Clayton, sclayton@uw.edu

% Riley's data
YD = [15,75,105,135,165, 255]; % yearday, assuming a simplified calendar with 12 30-day months
T = [4.61, 2.60, 3.81, 5.14, 9.66, 15.24]; % temperature
P = [33.7, 34.7, 21.9, 16.6, 19.2, 14.4]; % phosphate
N = [209, 172, 129, 285, 155, 153]; % nitrate
Z = [14, 24, 32, 106, 103, 135]; %  zooplankton, thousands per m^2
Ph = [118, 828, 2303, 871, 478, 560]; % phytoplankton pigments, thousands of Harvey units per m^2

% Make a plot of Riley's data
figure(1);
subplot(1,4,1)
plot(YD, Ph);ylabel('Phytoplankton');
subplot(1,4,2)
plot(YD, N);ylabel('Nitrate');
subplot(1,4,3)
plot(YD, T);ylabel('Temperature');
subplot(1,4,4)
plot(YD, Z);ylabel('Zooplankton');xlabel('Yearday')

% Riley's model
% Constant parameters
p = 2.5; % p, photosynthetic constant, 1/day
R0 = 0.0175; % R0, respiratory rate, 1/day
r = 0.069; % r, rate of change of respiratory rate with temperature, 1/C 
g = 0.0075; % g, grazing rate, 1/day 
P0 = 3.4; % P0, initial phytoplankton concentration, g C/m^2

% Seasonally varying parameters
% I0, incident solar radiation
I0 = [0.088, 0.094, 0.112, 0.138, 0.174, 0.212, 0.247, 0.272, 0.290, 0.306, 0.321, 0.329, 0.319, 0.302, 0.284, 0.267, 0.250, 0.230, 0.204, 0.174, 0.144, 0.115, 0.094, 0.086];
% k, extinction coefficient
k = [0.121, 0.121, 0.124, 0.128, 0.136, 0.145, 0.159, 0.2, 0.205, 0.17, 0.17, 0.17, 0.17, 0.17, 0.165, 0.162, 0.159, 0.154, 0.145, 0.138, 0.131, 0.126, 0.121, 0.121];
% z1, depth of euphotic zone
z1 = [34., 34., 35., 35., 35., 34., 32., 26., 25., 31., 32., 32., 32., 31., 32., 32., 32., 33., 34., 34., 35., 35., 34., 33.];
% 1-N, correction factor for nutrient depletion
N1 = [1.,1.,1.,1.,1.,1.,1.,1.,0.95,0.92,0.9,0.88,0.82,0.76,0.69,0.63,0.60,0.59,0.63,0.69,0.77,0.85,0.92,0.97];
% 1-V, correction factor for vertical turbulence
V1 = [0.64, 0.64, 0.69, 0.73, 0.78, 0.85, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.95, 0.76, 0.66];
% Ph, estimated mean photosynthetic rate
Ph1 = [0.034, 0.036, 0.044, 0.055, 0.071, 0.091, 0.12, 0.13, 0.131, 0.132, 0.134, 0.134, 0.122, 0.108, 0.093, 0.081, 0.073, 0.067, 0.065, 0.063, 0.06, 0.053, 0.039, 0.033];
% T, mean surface temperature
T = [5.2, 4.1, 3.2, 2.7, 2.4, 2.5, 2.7, 3.4, 4.5, 5.9, 7.6, 9.7, 11.8, 13.9, 15.5, 16.3, 16.6, 16.4, 15.5, 14.2, 12.4, 10.5, 8.5, 6.7];
% Rt, estimated phytoplankton respiration
Rt = [0.024, 0.023, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.024, 0.026, 0.03, 0.035, 0.038, 0.045, 0.051, 0.054, 0.056, 0.054, 0.051, 0.045, 0.04, 0.037, 0.031, 0.028];
% Z, mean zooplankton biomass
Z = [1.3, 1.4, 1.7, 2.6, 4.2, 5.8, 7.5, 8.9, 17.2, 19.3, 18.8, 14.0, 6.9, 6.2, 6.0, 5.7, 5.1, 4.5, 3.9, 3.3, 3.2, 2.6, 2.0, 1.6];
% G, grazing rate
G1 = [0.01, 0.011, 0.013, 0.02, 0.031, 0.043, 0.056, 0.067, 0.129, 0.145, 0.141, 0.105, 0.052, 0.047, 0.045, 0.043, 0.038, 0.034, 0.029, 0.025, 0.024, 0.02, .0015, 0.012];

% time stepping
dt = 1; % each time step is 1 day long
day = 1:360;

% interpolate the forcing data onto a daily timestep
I0d = interp1(day, 1:15:360, I0)
kd = interp1(day, 1:15:360, k)
z1d = interp1(day, 1:15:360, z1)
N1d = interp1(day, 1:15:360, N1)
V1d = interp1(day, 1:15:360, V1)
Td = interp1(day, 1:15:360, T)
Zd = interp1(day, 1:15:360, Z)

% pre-assign arrays for results of model
Pmod = zeros(360,1);
Phmod = zeros(360,1);
Rmod = zeros(360,1);
Gmod = zeros(360,1);

dP = zeros(360,1);
P(1) = P0;


% integrate the model forward in time
for j =1:359;
    Phmod(j) = (p*I0d(j))/(kd(j)*z1d(j)) * (1 - exp(-kd(j)*z1d(j))) * N1d(j) * V1d(j);
    Rmod(j) = R0 * exp(r*Td(j));
    Gmod(j) = g*Zd(j);
    dP(j) = Pmod(j)*(Phmod(j)-Rmod(j)-Gmod(j));
    Pmod(j+1) = P(j) + dP(j)*dt;

end

% plot the model result vs the observations

figure(2);
plot(day, Phmod, YD, P, 'or');ylabel('Phytoplankton');xlabel('Yearday')





