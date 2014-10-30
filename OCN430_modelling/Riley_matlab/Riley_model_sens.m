% Riley_model_sens script
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

% prompt user to pick which parameter to test
test = 0;
while test<1 || test>4,
test = Riley_input;
end

% create upper and lower bounds on the tested parameter
% these are +/- 20% of the original value

if test == 1
    P0_low = 0.8*P0;
    P0_up = 1.2*P0;
else
    P0_low=P0;
    P0_up=P0;
end

if test == 2
    p_low = 0.8*p;
    p_up = 1.2*p;
else
    p_low=p;
    p_up=p;
end

if test == 3
    R0_low = 0.8*R0;
    R0_up = 1.2*R0;
else
    R0_low=R0;
    R0_up=R0;
end

if test == 4
    g_low = 0.8*g;
    g_up = 1.2*g;
else
    g_low=g;
    g_up=g;
end

% pre-assign arrays for results of model
Pmod = zeros(360,1);
Phmod = zeros(360,1);
Rmod = zeros(360,1);
Gmod = zeros(360,1);

P_up = zeros(360,1);
Ph_up = zeros(360,1);
R_up = zeros(360,1);
G_up = zeros(360,1);

P_low = zeros(360,1);
Ph_low = zeros(360,1);
R_low = zeros(360,1);
G_low = zeros(360,1);

dP = zeros(360,1);
dP_up = zeros(360,1);
dP_low = zeros(360,1);

% set inital condition on phytoplankton
Pmod(1) = P0;
P_up(1) = P0_up;
P_low(1) = P0_low;

% integrate the model forward in time
for j =1:359;
    Phmod(j) = (p*I0d(j))/(kd(j)*z1d(j)) * (1 - exp(-kd(j)*z1d(j))) * N1d(j) * V1d(j);
    Rmod(j) = R0 * exp(r*Td(j));
    Gmod(j) = g*Zd(j);
    dP(j) = Pmod(j)*(Phmod(j)-Rmod(j)-Gmod(j));
    Pmod(j+1) = Pmod(j) + dP(j)*dt;
    
    Ph_up(j) = (p_up*I0d(j))/(kd(j)*z1d(j)) * (1 - exp(-kd(j)*z1d(j))) * N1d(j) * V1d(j);
    R_up(j) = R0_up * exp(r*Td(j));
    G_up(j) = g_up*Zd(j);
    dP_up(j) = P_up(j)*(Ph_up(j)-R_up(j)-G_up(j));
    P_up(j+1) = P_up(j) + dP_up(j)*dt;
    
    Ph_low(j) = (p_low*I0d(j))/(kd(j)*z1d(j)) * (1 - exp(-kd(j)*z1d(j))) * N1d(j) * V1d(j);
    R_low(j) = R0_low * exp(r*Td(j));
    G_low(j) = g_low*Zd(j);
    dP_low(j) = P_low(j)*(Ph_low(j)-R_low(j)-G_low(j));
    P_low(j+1) = P_low(j) + dP_low(j)*dt;

end
disp('Success!')

% plot the model result vs the observations
figure;
plot(day, Pmod,'-k', day, P_low ,'-.g',day, P_up, '-.b', YD, Ph*17/1000, 'or');
ylabel('Phytoplankton (g C m^{-2})');xlabel('Yearday')


