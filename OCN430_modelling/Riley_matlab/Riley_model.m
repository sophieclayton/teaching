% Riley_model script
% Sophie Clayton, 2014
% sclayton@uw.edu

clear all

load Riley_params
load Riley_data

% interpolate Riley forcing data onto daily timestep
iday = 0:15:360;
iday(1)=1;
day = 1:360;
dt = 1;

I0d = interp1(iday, I0, day);
kd = interp1(iday, k, day);
z1d = interp1(iday, z1, day);
N1d = interp1(iday, N1, day);
V1d = interp1(iday, V1, day);
Td = interp1(iday, T, day);
Zd = interp1(iday, Z, day);

% pre-assign arrays for results of model
Pmod = zeros(360,1);
Phmod = zeros(360,1);
Rmod = zeros(360,1);
Gmod = zeros(360,1);

dP = zeros(360,1);
Pmod(1) = P0;


% integrate the model forward in time
for j =1:359;
    Phmod(j) = (p*I0d(j))/(kd(j)*z1d(j)) * (1 - exp(-kd(j)*z1d(j))) * N1d(j) * V1d(j);
    Rmod(j) = R0 * exp(r*Td(j));
    Gmod(j) = g*Zd(j);
    dP(j) = Pmod(j)*(Phmod(j)-Rmod(j)-Gmod(j));
    Pmod(j+1) = Pmod(j) + dP(j)*dt;

end
disp('Success!')
% plot the model result vs the observations

figure(2);
plot(day, Pmod, YD, Ph*17/1000, 'or');ylabel('Phytoplankton (g C m^{-2})');xlabel('Yearday')


