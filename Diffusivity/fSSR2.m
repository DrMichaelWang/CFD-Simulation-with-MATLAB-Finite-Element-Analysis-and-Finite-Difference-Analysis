function f = fSSR2(a,xm,ym)
sums = 0;
R = 8.31434;
T = 298.15;
VIL = 25;
rhoIL = 0.9;
MWIL = 18;
Vg = 30-VIL;
L = 4;
k = (8*R*T*VIL*rhoIL)/(pi^2*Vg*MWIL);
for n=0:55
    sums = sums + (1/(2*n+1)^2*exp((2*n+1)^2*pi^2*a(2)*xm)/(4*L^2)-1);
end
yp = (k/a(1))*sums;
f = sum((ym-yp).^2);