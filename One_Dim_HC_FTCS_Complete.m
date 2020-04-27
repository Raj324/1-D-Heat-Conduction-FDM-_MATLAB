L=input('Please specify the Length of ROD in mm.: ');
dx=input('Please specify the Space increments required along the Rod.: ');
tt=input('Please specify the No. of iterations required in Time Level: ');
al=input('Please specify the Thermal diffusivity of the Material of the Rod in mm2/sec: ');
T0X0=input('Please specify the Temperature of the Rod @ x=0 and t=0:');
T0XL=input('Please specify the Temperature of the Rod @ x=L and t=0:');
TX0=input('Please specify the Temperature of the Rod @ x=0 and t > 0:');
TXL=input('Please specify the Temperature of the Rod @ x=L and t > 0:');
N=(L/dx)+1;%Number of Grid Points%
dt=1;%Increments in time level set as 1%
x=0:dx:L;
vc=(al*dt)/(dx^2); %Setting up the value of Constant used in FTCS Analysis%
if vc<=0.5 %CHECKING THE STABILITY CRITERION%
    fprintf('STABLE- It is safe to continue \n');
for m=1:tt %Setting up the array for time level as per the input for No. of iterations%
    M(1,1)=0;
    M(m,1)=M(1,1)+((m-1)*dt);
end
for k=1:N %Setting up the array for space level as per the input for No. space increments along the rod%
    X(1,1)=0;
    X(1,k)=X(1,1) + ((k-1)*dx);
end
for i=2:N-1
    T(1,1)=T0X0;%Initializing the boundary condition at x=0 and at time,t= 0%
    T(1,N)=T0XL;%Initializing the boundary condition at x=L and at time,t= 0%
    T(1,i)=T(1,1)-(x(i)/L)*(T(1,1)-T(1,N));%temperature distribution along the Rod at time t=0, using Fouriers law of heat conduction%
end
for j=2:tt
    T(j,1)=TX0;%Initializing the boundary condition at x=0 and at time > 0%
    T(j,N)=TXL;%Initializing the boundary condition at x=L and at time > 0%
end
for j=2:tt %Using FTCS Method to solve for the Temperature distribution along the Rod using the BCs.%
    i=2:N-1;
    T(j,i)=vc*(T(j-1,i-1) + T(j-1,i+1)-(2*T(j-1,i))) + T(j-1,i);
end
B=[zeros(1,1);M];
C=[X;T];
A=[B,C]; % Setting up Array to display Space Increment in along the First row except in the (1,1) cell, Time level along the First Column, and The temp Distribution of the corresponding space level and the respective time level%
J=[T(1,:);T(11,:);T(21,:);T(31,:);T(41,:);T(51,:)];
Mid_Temp=T(:,(N+1)/2); % Array for Mid Point Temperature at all time levels %
    disp(A)
    plot(X,T)%Plot of x v/s T%
    title('Variations in Temperature along the Rod in Different time levels')
    xlabel('Space Increments along the Rod in mm') % x-axis label
    ylabel('Temperature in Deg Celcius') % y-axis label
    figure
    plot(X,J)%Plot of x v/s T for time level t=0,10,20,30,40,50 sec%
    title('Variations in Temperature along the Rod in time level t=1,10,20,30,40,50 Sec')
    xlabel('Space Increments along the Rod in mm') % x-axis label
    ylabel('Temperature in Deg Celcius') % y-axis label
    figure
    plot(M,Mid_Temp) % Plot of t vs Mid Point Temp (T[x/2]) %
    xlabel('Time Level in seconds') % x-axis label
    ylabel('Rod Mid Point Temperature in Deg Celcius') % y-axis label
else
    fprintf('Warning!!!! - Unstability detected and your PC may blow up ');
end
xlswrite('D:\CFD ASSIGN1\TempFTCS.xlsx',A) % Writing the calculations in an Excel file%