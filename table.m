%%
clear
%%
vinf=2.5;
[ene_res, rp_res, ra_res]=TG_resonant('Ganimede', vinf);

for i=1:10
    for j=1:10
        if rp_res(j,i) ~= 0
            [alpha(j,i), ~]=TG_get_alpha_ene('Ganimede', rp_res(j,i), vinf);
        end
    end
end
%%
for i = 1:10
    for j=1:10
        alpha_v(10*(i-1) + j) = real(alpha(i,j));
    end
end
delta_max = get_deltamax('Ganimede', vinf);
delta_max = delta_max *pi/180;
for i = 1:100
    for j=1:100
        delta(i,j) = abs(alpha_v(i) - alpha_v(j));
        if delta(i,j) > delta_max
            delta(i,j) = -100;
        end
    end
end
%%
for i = 1:100
    for j=1:100
        if delta(i,j) > 0
            lon(i,j) = ((alpha_v(i) + pi + alpha_v(j)) / 2) * 180/pi - 270;
        end
    end
end

%%% PROBLEMA ALPHA MAI > PI %%%

%% Sections
p = (1+sqrt(5))/2;
b1 = [3*p, 2*p, -p, -2*p, -3*p];
b2 = [-1, -(2+p), -(1+2*p), -(2+p), -1];
sect = [8,30,22,21,12,5];
for i = 1:length(b1)
    lim(i+1) = atan2(b2(i), b1(i)) * 180/pi;
end
lim(length(b1)+2)=-180;
l = length(lim);
l = l-1;
t=[];
for i = 1:100
    for j=1:100
        for k = 1:l
            if lon(i,j) <= lim(k) && lon(i,j) > lim(k+1) && lon(i,j)~=0
                T(i,j) = sect(k);
                t(k) = sect(k);
            end
        end
    end
end

        