function [alpha, v_inf]=TG_get_alpha_vinf(moon,peri,apo)
% The function computes the values of alpha (pump angle) and v_inf (hyperbolic eccess velocity) of the s/c flying on a given orbit (periapsis and apoapsis are inputs) around a desired moon
% Also, it plots the Tisserand curve of interest highlighting the input orbit

moons=["Io", "Europa", "Ganimede", "Callisto"];
Rp = [421.6e3 670.9e3 1.07e6 1.883e6]; % moon distance -> Io, Europa, Ganimede and Callisto respectevely
mu_jup = 1.899*10^27 * 6.6743 * 10^(-20); %Jupiter gravitational constant
for i=1:4
    if moons(i) == moon
        Rm = Rp(i);
    end
end
Vp = sqrt(mu_jup/Rm);   % Moon's velocity around Jupiter

% computing the orbital parameters of the orbit
a = (apo+peri)/2;
ene = -mu_jup /(2*a);
e = 1-peri/a;

% Computing the v_inf from a and e
v_inf = Vp *sqrt(3-Rm/a- sqrt((1-e^2)*4*a/Rm));


% Plotting the Tisserand curve
alpha1 = linspace(0,pi);
v1 = sqrt(Vp^2 + v_inf^2 + 2*Vp*v_inf.*cos(alpha1));
a1 = -mu_jup/2./(v1.^2/2 - mu_jup/Rm);
e1 = sqrt(1-Rm./a1.*(0.5*(3-Rm./a1-(v_inf./Vp).^2)).^2);
rp=a1.*(1-e1);
ra=a1.*(1+e1);
ene1=-mu_jup./(2.*a1);

figure(1)
set(gca, 'XScale','log', 'YScale', 'log')
hold on
grid on
plot(rp,ra)
xlabel('rp')
ylabel('ra')
figure(2)
set(gca, 'XScale','log')
hold on
grid on
plot(rp,ene1)
xlabel('rp')
ylabel('orbital energy')

% Computing alpha from the geometric relation of the velocities
v1 = sqrt(mu_jup*2/Rm - mu_jup/a);
alpha = acos((v1^2 - Vp ^2 -v_inf^2)/(2*Vp*v_inf));


% Highlighting the input orbit
figure(1)
hold on
plot(peri,apo,'s','Color','r')
figure(2)
hold on
plot(peri, ene, 's', 'Color','r')
