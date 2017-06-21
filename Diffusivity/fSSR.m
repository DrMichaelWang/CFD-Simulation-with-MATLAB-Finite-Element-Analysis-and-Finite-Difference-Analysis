function f = fSSR(a,xm,ym)
sums = 0;
R = 8.31434;
T = input('Enter Temperature (K):');
VIL = input('Enter Lubricant Volume (cm3):');
rhoIL = input('Enter Lubricant Density (g/cm3):');
MWIL = input('Enter Lubricant Molecular Weight (g/mol):');
Vg = 30-VIL;
L = input('Enter Lubricant Depth (cm):');
k = (8*R*T*VIL*rhoIL)/(pi^2*Vg*MWIL);
for n=0:55
    sums = sums + (1/(2*n+1)^2*exp((2*n+1)^2*pi^2*a(2)*xm)/(4*L^2)-1);
end
yp = (k/a(1))*sums;
f = sum((ym-yp).^2);