% JOVIAN CAPTURE
clear
Tg_JovianMoons
clear
v_first = 3.4; % km/s
R_jup = 69911; % km  Jupiter radius
r_first = 1000 * R_jup; % km
mu_jup = 1.899*10^27 * 6.6743 * 10^(-20); %Jupiter gravitational constant

ene = v_first^2 / 2 - mu_jup / r_first;

peri = 1.1e6;
a = - mu_jup / (2*ene);
apo = 2*a-peri;

% CALLISTO FLYBY
% Pre-flyby
[alpha0, vinf] = TG_get_alpha_vinf('Callisto', peri, apo); 
TG_resonant('Callisto', vinf);
peri1=0.98e6;
[alpha1, ene]=TG_get_alpha_ene('Callisto',peri1,vinf);
delta_max = get_deltamax('Callisto',vinf);
% Compute alpha diff
if abs((alpha1-alpha0)*180/pi) <= delta_max
    disp('Flyby possible')
else
    disp('Fliby not achievable')
end
% Post-flyby
a1 = -mu_jup/(2*ene);
apo1 = 2*a1-peri1;
[alpha2, vinf2] = TG_get_alpha_vinf('Ganimede', peri1, apo1);
TG_resonant('Ganimede', vinf2)
