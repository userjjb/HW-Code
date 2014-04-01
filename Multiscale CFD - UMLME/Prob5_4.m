clear all
close all

choice =4; %Choose part a/b and flux scheme

switch choice
case 1%Part a (central Difference)-----------------------------------------
    phi_initial = 1;
    phi_final = 0;

    K = 100;                %Number of cells
    x = 0:1/(K+1):1;        %Cell center locations
    xi = 1/(2*(K+1)):1/(K+1):1; %Cell interfaces
    deltax = x(2)-x(1);     %Equidistant grid, constant deltax

    mxi = 40*(1-xi);        %Mass flow rate at interfaces

    A = zeros(K+2);
    for i=2:length(A)-1     %Assemble linear eqn matrix
        A(i,i-1:i+1) = [-(1/deltax + (mxi(i-1)/2)) ...
                        (2/deltax)+(mxi(i)-mxi(i-1))/2 ...
                        -(1/deltax - (mxi(i)/2))];
    end

    %Apply forcing
    b = zeros(length(A),1);

    %Apply BCs
    A(1,:) = 0; A(1,1) = 1;
    b(1)= phi_initial;
    A(end,:) = 0; A(end,end) = 1;
    b(end) = phi_final;

    u=A\b;

    plot(x,u)
case 2%Part a (Power Law)--------------------------------------------------
    phi_initial = 1;
    phi_final = 2;

    K = 100;                %Number of cells
    x = 0:1/(K+1):1;        %Cell center locations
    xi = 1/(2*(K+1)):1/(K+1):1; %Cell interfaces
    deltax = x(2)-x(1);     %Equidistant grid, constant deltax

    mxi = 40*(1-xi);        %Mass flow rate at interfaces

    A = zeros(K+2);
    for i=2:length(A)-1     %Assemble linear eqn matrix
        A(i,i-1:i+1) = [-(((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))+mxi(i-1)) ...
                        (((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))+mxi(i-1))+((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))+ (mxi(i)-mxi(i-1)) ...
                        -((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))];
    end

    %Apply forcing
    b = zeros(length(A),1);

    %Apply BCs
    A(1,:) = 0; A(1,1) = 1;
    b(1)= phi_initial;
    A(end,:) = 0; A(end,end) = 1;
    b(end) = phi_final;

    u=A\b;

    plot(x,u)
case 3%Part b (central Difference)-----------------------------------------
    phi_initial = 1;
    phi_final = 0;
    mL = -40;

    K = 100;                %Number of cells
    x = 0:1/(K+1):1;        %Cell center locations
    xi = 1/(2*(K+1)):1/(K+1):1; %Cell interfaces
    deltax = x(2)-x(1);     %Equidistant grid, constant deltax

    mxi = 40*(1-xi);        %Mass flow rate at interfaces

    A = zeros(K+2);
    for i=2:length(A)-1     %Assemble linear eqn matrix
        A(i,i-1:i+1) = [-(1/deltax + (mxi(i-1)/2)) ...
                        (2/deltax)+(mxi(i)-mxi(i-1))/2 ...
                        -(1/deltax - (mxi(i)/2))];
    end
    
    for i=2:length(A)-1     %Apply linear term
        A(i,i) = A(i,i) + -mL;
        A(1,i) = A(1,i) - -mL;
    end

    %Apply forcing
    b = zeros(length(A),1);

    %Apply BCs
    A(1,:) = 0; A(1,1) = 1;
    b(1)= phi_initial;
    A(end,:) = 0; A(end,end) = 1;
    b(end) = phi_final;

    u=A\b;

    plot(x,u)
case 4%%Part b (Power Law)-------------------------------------------------
    phi_initial = 1;
    phi_final = 0;
    mL = -40;

    K = 100;                %Number of cells
    x = 0:1/(K+1):1;        %Cell center locations
    xi = 1/(2*(K+1)):1/(K+1):1; %Cell interfaces
    deltax = x(2)-x(1);     %Equidistant grid, constant deltax

    mxi = 40*(1-xi);        %Mass flow rate at interfaces

    A = zeros(K+2);
    for i=2:length(A)-1     %Assemble linear eqn matrix
        A(i,i-1:i+1) = [-(((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))+mxi(i-1)) ...
                        (((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))+mxi(i-1))+((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))+ (mxi(i)-mxi(i-1)) ...
                        -((1/deltax)*max(0,(1-(mxi(i-1)/(10*deltax)))^5))];
    end

    for i=2:length(A)-1     %Apply linear term
        A(i,i) = A(i,i) + -mL;
        A(1,i) = A(1,i) - -mL;
    end
    %Apply forcing
    b = zeros(length(A),1);

    %Apply BCs
    A(1,:) = 0; A(1,1) = 1;
    b(1)= phi_initial;
    A(end,:) = 0; A(end,end) = 1;
    b(end) = phi_final;

    u=A\b;

    plot(x,u)
end