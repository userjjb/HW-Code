%Josh Bevan 2014, 22.559 MS-CFD
%Solves a 1D steady-state conduction problem with periodic boundary
%conditions with heat generation 

%________Input Params_______________________
Bi= 0.05;   %Biot number for system
L= 20;      %Length of rod in diameters
xH= 3.5;    %Location of start of heated section, use n+1/2 to line up CV
LH= 2;      %Length of heated section
CV= 20;     %Number of CVs to use
d= 1/(L/CV)^2;%Grid size

%The primarily tridiagonal structure is a simple Toeplitz matrix; Including
%the same node twice in the domain necessitates a modification so the 2nd
%derivative jumps the coincident node
Theta = ( toeplitz([-2*d-4*Bi, d, 2:CV<0])...
        +fliplr(blkdiag([0,d],zeros(CV-1,CV-3),[d,0])) )...
        \ -and(0:L/CV:L>xH,0:L/CV:L<xH+LH)'


plot(0:L/CV:L,Theta)
xlabel('Axial Location (x/D)')
ylabel('Non-dimensional Temperature (theta)')