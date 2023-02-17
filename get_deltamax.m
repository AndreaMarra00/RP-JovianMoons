function delta_max = get_deltamax(moon, v_inf)
moons=["Io", "Europa", "Ganimede", "Callisto"];
Rm=[3642/2 3120/2 5268/2 4800/2];
mu = [8.932e22 4.791e22 1.482e23 1.077e23] * 6.6743 * 10^(-20);
h_min = 50;
for i=1:4
    if moons(i)==moon
        Rm = Rm(i);
        mu = mu(i);
    end
end

delta_max=2*asin((1+(Rm+h_min)*v_inf^2/mu)^(-1))*180/pi;