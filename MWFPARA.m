function v=MWFPARA

rou=980; %rou is the density of CO2,which unit is kg/m3
mu=10^(-1);%mu is the dynamic viscosity of soybean oil,which unit is Pa*s (kg/(m*s))
Vy0=5;%Vy0 is the initial velocity, which unit is m/s
Vc=-45.7/60; %Vc is the cutting speed, which unit is m/s
g=9.81; %g is the standard gravity, which unit is m/s^2
L1=0.4149*10^(-3); %Dimension along the X direction,its unit is m
L2=4.74*10^(-3); %Dimension along the Y direction,its unit is m
a=L1; %a is the constant term in the L(y) function
b=-0.0875; %b is the constant term in the L(y) function
xmax=200;
dx=L1/xmax;
dy=dx;
ymax=round(L2/dy)

%Initialize the boundary conditions
for i=2:xmax
    v(i,1)=Vy0;
    v(i,ymax+1)=0;
end

for j=1:ymax+1
    v(1,j)=Vc;
    v(xmax+1,j)=0;
end
%Give Initial Values
for i=2:xmax
    for j=2:ymax
        v(i,j)=0.99*v(i,j-1);
    end
end

%Iterations for Calculation Velocity Profile
for k=1:40000
    xstartpoint=xmax+1;
  for j=2:ymax
       for i=2:xstartpoint-1
        TERM1=(mu/rou)*(v(i,j+1)+v(i,j-1)+v(i+1,j)+v(i-1,j))/(dy^2);
        TERM2=-1*g;
        TERM3=0;
        TERM4=(v(i,j+1)-v(i,j-1))/(2*dy);
        TERM5=mu*4/(rou*(dy^2));
        v(i,j)=(TERM1+TERM2+TERM3)/(TERM4+TERM5);
        if v(i,j)<0, v(i,j)=0; end
       end
  end
end
