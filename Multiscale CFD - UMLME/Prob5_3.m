Hin=100;        %Input temps
Cin=20;         

N = 10;         %Number of cells
X = 0:1/N:1;    %Cell center locations
dX=1/N;         %Spatial step size

H = zeros(length(X),1);     %Cell averages
H(1)= Hin;                  %Input temps
C = zeros(length(X),1);
C(1)= Cin;

dH = zeros(length(X),1);    %Temp gradient

for i=1:length(H)-1
    dH = C(i)-H(i);         %Based on DE
    H(i+1) = H(i)+ dH*dX;   %Using upwind flux for Peclet=infinity
    C(i+1) = C(i)- dH*dX;   %Conservation of energy
end

uH=(1+exp(-2.*X))/2;         %Analytical solutions
uC=(1-exp(-2.*X))/2;

plot(X,(H-Cin)/(Hin-Cin),'--r')
hold on
plot(X, uH,'r')
plot(X,(C-Cin)/(Hin-Cin),'--b')
plot(X, uC,'b')
xlabel('x/L')
ylabel('Dimensionless Temperature')
title(['Parallel Flow Heat Exchanger, n=' num2str(N)''])
legend('Numerical','Analytical')