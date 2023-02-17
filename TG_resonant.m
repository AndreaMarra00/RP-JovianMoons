function [ene_list, rp_list, ra_list] = TG_resonant(moon, v_inf)
% The function plots and lists the main resonant orbits with the desired moon given the v_inf
% The function's outputs are three 3x3 matrixes:
    % rp_list : a matrix containing the periapsis' values for the various resonant orbits found: rp_list(j,k) corresponds to the periapsis of the orbit with a resonance of j:k
    % ra_list : a matrix containing the apoapsis' values for the various resonant orbits found: rp_list(j,k) corresponds to the apoapsis of the orbit with a resonance of j:k
    % ene_list : a matrix containing the orbital energy values for the various resonant orbits found: ene_list(j,k) corresponds to the orbital energy of the orbit with a resonance of j:k
%provaprovaprova 
provaprovaprova
moons=["Io", "Europa", "Ganimede", "Callisto"];
Rp = [421.6e3 670.9e3 1.07e6 1.883e6]; % moon distance -> Io, Europa, Ganimede and Callisto respectevely
mu_jup = 1.899*10^27 * 6.6743 * 10^(-20); %Jupiter gravitational constant
for i=1:4
    if moons(i) == moon
        Rm = Rp(i);
    end
end
Vp = sqrt(mu_jup/Rm);   % Moon's velocity around Jupiter
Tp = 2*pi * sqrt(Rm^3 / mu_jup);    % Moon's orbital period


% Plotting the Tisserand curve linked to the inputs
alpha=linspace(0,pi);
v1 = sqrt(Vp^2 + v_inf^2 + 2*Vp*v_inf.*cos(alpha));
a = -mu_jup/2./(v1.^2/2 - mu_jup/Rm);
e = sqrt(1-Rm./a.*(0.5*(3-Rm./a-(v_inf./Vp).^2)).^2);
rp=a.*(1-e);
ra=a.*(1+e);
ene=-mu_jup./(2.*a);
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
plot(rp,ene)
xlabel('rp')
ylabel('orbital energy')

% Plotting the resonances
for k = 1:3 % multiple of the moons orbital period
    for j = 1:3 % multiple of the s/c orbital period
        a_res = (k * Tp * sqrt(mu_jup) / ( 2 * pi * j)) ^ (2/3); % Semimajor axis in function of multiples of the moons period 
        e_res = sqrt(1-Rm/a_res*(0.5*(3-Rm/a_res-(v_inf/Vp)^2))^2);
        alpha_res(k,j)=acos((Vp^2*(2-(j/k)^(2/3))-v_inf^2-Vp^2)/(2*v_inf*Vp));

        ene_res=-mu_jup/(2*a_res);
        rp_res = a_res*(1-e_res);
        ra_res = a_res*(1+e_res);

        if imag(alpha_res(k,j)) == 0
            % Creating the two output matrixes 
            rp_list(j,k) = rp_res;
            ra_list(j,k) = ra_res;
            ene_list(j,k) = ene_res;

            figure(1)
            hold on
            %                 M=["" "o" "V"; "d" "" "s"; "^" "<" ""];%Matrix of the markers to see the resonance points
            m(k,j)=plot(rp_res, ra_res,'o','MarkerSize',4);  
            figure(2)
            hold on
            m1(k,j)=plot(rp_res, ene_res,'o','MarkerSize',4);
        end
    end
end 

