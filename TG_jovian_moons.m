clear, close all
%% TG of Jovian Moons

Rp = [421.6e3 670.9e3 1.07e6 1.883e6]; % moon distance -> Io, Europa, Ganimede and Callisto respectevely
Rm=[3642/2 3120/2 5268/2 4800/2]; %radius of the moons: Io, Europa, Ganimede and Callisto in km
colors = ['r' 'b' 'g' 'c'];
h_min=50; %minimum flyby altitude in km
mu_jup = 1.899*10^27 * 6.6743 * 10^(-20); %Jupiter gravitational constant
moons=["Io", "Europa", "Ganimede", "Callisto"];
Tp = 2*pi * sqrt(Rp.^3 / mu_jup);  % Moons orbital period
M=["p" "o" ""; "d" "p" "s"; "^" "<" "p"];%Matrix of the markers to see the resonance points

% moon speed km/s -> Io, Europa, Ganimede and Callisto respectevely:
Vp = [];  
for j=1:4
    Vp(j) = sqrt(mu_jup/Rp(j));
end
%% PERIAPSIS vs APOAPSIS
alpha=linspace(0,pi);
for i=1:4
    v_inf = 1;  %Initialized v infinite km/s (example)
    while v_inf <= 7
        v1 = sqrt(Vp(i)^2 + v_inf^2 + 2*Vp(i)*v_inf.*cos(alpha));
        a = -mu_jup/2./(v1.^2/2 - mu_jup/Rp(i));
        e = sqrt(1-Rp(i)./a.*(0.5*(3-Rp(i)./a-(v_inf./Vp(i)).^2)).^2);
        rp=a.*(1-e);
        ra=a.*(1+e);
        ene=-mu_jup./(2.*a);
        figure(1)
        set(gca, 'XScale','log', 'YScale', 'log')
        hold on
        grid on
        if i==1
            p1=plot(rp,ra,colors(i));
        elseif i==2
            p2=plot(rp,ra,colors(i));
        elseif i==3
            p3=plot(rp,ra,colors(i));
        elseif i==4
            p4=plot(rp,ra,colors(i));
        end
        xlabel('R_p')
        ylabel('R_a')
        title('Apoapsis vs Periapsis')
        m=[];   
        for k = 1:3 % multiple of the moons orbital period
            for j = 1:3 % multiple of the s/c orbital period
                a_res = (k * Tp(i) * sqrt(mu_jup) / ( 2 * pi * j)) ^ (2/3); % Semimajor axis in function of multiples of the moons period 
                e_res = sqrt(1-Rp(i)/a_res*(0.5*(3-Rp(i)/a_res-(v_inf/Vp(i))^2))^2);
                alpha_res(k,j)=acos((Vp(i)^2*(2-(j/k)^(2/3))-v_inf^2-Vp(i)^2)/(2*v_inf*Vp(i)));
                rp_res = a_res*(1-e_res);
                ra_res = a_res*(1+e_res);
                if imag(alpha_res(k,j)) == 0
                    figure(1)
                    hold on
                    m(k,j)=plot(rp_res, ra_res, M(k,j),'MarkerSize',4,'MarkerFaceColor',colors(i),'MarkerEdgeColor',colors(i));  
                end
            end
        end 
        v_inf = v_inf + 2;
    end
end
legend([p1,p2,p3,p4,m(1,1),m(1,2),m(2,1),m(2,2),m(2,3),m(3,1),m(3,2),m(3,3)],'Io', 'Europa', 'Ganimede','Callisto','1:1','1:2','2:1','2:2','2:3','3:1','3:2','3:3')

%% Orbital energy vs periapsis
alpha=linspace(0,pi);
for i=1:4
    v_inf = 1;  %Initialized v infinite km/s (example)
    while v_inf <= 7
        v1 = sqrt(Vp(i)^2 + v_inf^2 + 2*Vp(i)*v_inf.*cos(alpha));
        a = -mu_jup/2./(v1.^2/2 - mu_jup/Rp(i));
        e = sqrt(1-Rp(i)./a.*(0.5*(3-Rp(i)./a-(v_inf./Vp(i)).^2)).^2);
        rp=a.*(1-e);
        ene=-mu_jup./(2.*a);
        figure(2)
        set(gca, 'XScale','log')
        hold on
        grid on
        if i==1
            e1=plot(rp,ene,colors(i));
        elseif i==2
            e2=plot(rp,ene,colors(i));
        elseif i==3
            e3=plot(rp,ene,colors(i));
        elseif i==4
            e4=plot(rp,ene,colors(i));
        end
        xlabel('R_p')
        ylabel('Orbital Energy')
        title('Orbital Energy vs Periapsis')
        m_ene=[];   
        for k = 1:3 % multiple of the moons orbital period
            for j = 1:3 % multiple of the s/c orbital period
                a_res = (k * Tp(i) * sqrt(mu_jup) / ( 2 * pi * j)) ^ (2/3); % Semimajor axis in function of multiples of the moons period 
                e_res = sqrt(1-Rp(i)/a_res*(0.5*(3-Rp(i)/a_res-(v_inf/Vp(i))^2))^2);
                alpha_res(k,j)=acos((Vp(i)^2*(2-(j/k)^(2/3))-v_inf^2-Vp(i)^2)/(2*v_inf*Vp(i)));
                ene_res=-mu_jup/(2*a_res);
                rp_res = a_res*(1-e_res);
                if imag(alpha_res(k,j)) == 0
                    figure(2)
                    hold on
                    m_ene(k,j)=plot(rp_res, ene_res, M(k,j),'MarkerSize',4,'MarkerFaceColor',colors(i),'MarkerEdgeColor',colors(i));
                end
            end
        end 
        v_inf = v_inf + 2;
    end
end
legend([e1,e2,e3,e4,m_ene(1,1),m_ene(1,2),m_ene(2,1),m_ene(2,2),m_ene(2,3),m_ene(3,1),m_ene(3,2),m_ene(3,3)],'Io', 'Europa', 'Ganimede','Callisto','1:1','1:2','2:1','2:2','2:3','3:1','3:2','3:3')


%% Maximum Deflection angle
delta_max=[];
mu = [8.932e22 4.791e22 1.482e23 1.077e23] * 6.6743 * 10^(-20);
for i = 1:4
    j=1;
    for v_inf=1:2:7 
        delta_max(i,j)=2*asin((1+(Rm(i)+h_min)*v_inf^2/mu(i))^(-1))*180/pi;
        j=j+1;
    end
end

% print and save in the directory the maximum deflection angle table 
fid=fopen('Table_delta.txt','w');
fprintf(fid,'Maximum Deflection Angle(°) vs V_infinity(km/s) table:\n');
fprintf(fid,'Moons\t |\t V_inf=1\t |\t V_inf=3\t |\t V_inf=5\t |\t V_inf=7\n');
fprintf(fid,'Io\t         %6.3f\t         %6.3f\t         %6.3f\t        %6.3f\n', [delta_max(1,1);delta_max(1,2);delta_max(1,3);delta_max(1,4)]);
fprintf(fid,'Europa\t     %6.3f\t         %6.3f\t        %6.3f\t        %6.3f\n', [delta_max(2,1);delta_max(2,2);delta_max(2,3);delta_max(2,4)]);
fprintf(fid,'Ganimede\t %6.3f\t     %6.3f\t         %6.3f\t        %6.3f\n', [delta_max(3,1);delta_max(3,2);delta_max(3,3);delta_max(3,4)]);
fprintf(fid,'Callisto\t %6.3f\t         %6.3f\t         %6.3f\t        %6.3f\n', [delta_max(4,1);delta_max(4,2);delta_max(4,3);delta_max(4,4)]);
fclose(fid);





