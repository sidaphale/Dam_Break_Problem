%Godunov scheme for 1st order accurate
clear all; close all; clc;
N=101; 
L=100; 
delta_x=L/(N-1); 
g=9.81; 
delta_t=0.001; 
t=4; 
c=delta_t/delta_x; 
% intial values of h and v
for i=1:N
    for l=1
        if i*delta_x<50
            Qold(1,i)=10;
            Qold(1+1,i)=0;
        else
            Qold(1,i)=1;
            Qold(1+1,i)=0;
        end
    end
end
tinitial=0;
while tinitial<=t
    tinitial=tinitial+delta_t;
    for i=1:N
        c(i)=sqrt(g.*Qold(1,i));
    end
    %Calculating eigen values
    for i=1:N
        eigen(1,i)=(Qold(2,i)./Qold(1,i))+c(i);
        eigen(2,i)=(Qold(2,i)./Qold(1,i))-c(i);
    end
    alpha=max(abs(eigen));
    for i=1:N
        Eold(1,i)=Qold(2,i);
        Eold(2,i)=((Qold(2,i)^2/Qold(1,i))+(0.5*g*(Qold(1,i)^2)));
    end
    for i=1:N-1
        F1(1,i)=0.5*(Eold(1,i)+Eold(1,i+1))-0.5*max(alpha(i),alpha(i+1))*(Qold(1,i+1)-Qold(1,i));
        F1(2,i)=0.5*(Eold(2,i)+Eold(2,i+1))-0.5*max(alpha(i),alpha(i+1))*(Qold(2,i+1)-Qold(2,i));
    end
    for i=2:N
        F2(1,i)=0.5*(Eold(1,i)+Eold(1,i-1))-0.5*max(alpha(i),alpha(i-1))*(Qold(1,i)-Qold(1,i-1));
        F2(2,i)=0.5*(Eold(2,i)+Eold(2,i-1))-0.5*max(alpha(i),alpha(i-1))*(Qold(2,i)-Qold(2,i-1));
    end
    for i=2:N-1
        Qnew(1,i)=Qold(1,i)-(delta_t/delta_x)*(F1(1,i)-F2(1,i));
        Qnew(2,i)=Qold(2,i)-(delta_t/delta_x)*(F1(2,i)-F2(2,i));
    end
    for l=1:2
        i=1;
        Qnew(l,i)=Qold(l,i)-(delta_t/delta_x)*(F1(l,i+1)-F1(l,i));
    end
    for l=1:2
        i=N;
        Qnew(l,i)=Qold(l,i)-(delta_t/delta_x)*(F2(l,i)-F2(l,i-1));
    end
    Qold=Qnew;
    for i=1:N
        Unew(1,i)=Qnew(2,i)/Qold(1,i);
    end
figure(1);
plot(Qnew(1,:));
title('Plot for flow of water using Godunov scheme');
xlabel('x in m'); ylabel('H in m');
pause(0.0001);
end
figure(1);
plot(Qnew(1,:));
title('Plot for flow of water using Godunov scheme');
xlabel('x in m'); ylabel('H in m');
figure(2);
plot(Unew(1,:));
title('Plot for velocity of water');
xlabel('x in m'); ylabel('U in m/s');        