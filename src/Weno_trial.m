%WENO Scheme
function Weno_trial
clear all; close all; clc;
N=101; 
L=100; 
delta_x=L/(N-1);
g=9.81; 
delta_t=0.001; 
t=4; 

% intial values of h and v
for i=1:N
    for l=1
        if i*delta_x<50
            Qold(l,i)=10;
            Qold(l+1,i)=0;
        else
            Qold(l,i)=1;
            Qold(l+1,i)=0;
        end
    end
end
tinitial=0;
z=0;
while tinitial<=t
    tinitial=tinitial+delta_t;
    z=z+1;
    [delta_Eold]=Weno(Qold, delta_t);
    for l=1:2
        for i=4:N-2
            Q2(l,i)=Qold(l,i)-(1/4)*(delta_t/delta_x)*delta_Eold(l,i);
        end
        Q2(l,1:3)=Q2(l,4);
        Q2(l,N-1:N)=Q2(l,N-2);
    end
    for i=1:N
        c(i)=sqrt(g.*Q2(1,i));
    end
    %Calculating eigen values
    for i=1:N
        eigen(1,i)=(Q2(2,i)./Q2(1,i))+c(i);
        eigen(2,i)=(Q2(2,i)./Q2(1,i))-c(i);
    end
    d(z)=max(max(abs(eigen)));
    [delta_E1]=Weno(Q2,delta_t);
    for l=1:2
        for i=4:N-2
            Q3(l,i)=Qold(l,i)-(1/3)*(delta_t/delta_x)*delta_E1(l,i);
        end
        Q3(l,1:3)=Q3(l,4);
        Q3(l,N-1:N)=Q3(l,N-2);
    end
    for i=1:N
        c(i)=sqrt(g.*Q3(1,i));
    end
    for i=1:N
        eigen(1,i)=(Q3(2,i)./Q3(1,i))+c(i);
        eigen(2,i)=(Q3(2,i)./Q3(1,i))-c(i);
    end
    e(z)=max(max(abs(eigen)));
    [delta_E2]=Weno(Q3,delta_t);
    for l=1:2
        for i=4:N-2
            Q4(l,i)=Qold(l,i)-(1/2)*(delta_t/delta_x)*delta_E2(l,i);
        end
        Q4(l,1:3)=Q4(l,4);
        Q4(l,N-1:N)=Q4(l,N-2);
    end
    for i=1:N
        c(i)=sqrt(g.*Q4(l,i));
    end
    for i=1:N
        eigen(1,i)=(Q4(2,i)./Q4(1,i))+c(i);
        eigen(2,i)=(Q4(2,i)./Q4(1,i))-c(i);
    end
    f(z)=max(max(abs(eigen)));
    [delta_E3]=Weno(Q4,delta_t);
    for l=1:2
        for i=4:N-2
            Qnew(l,i)=Qold(l,i)-(delta_t/delta_x)*delta_E3(l,i);
        end
        Qnew(l,1:3)=Qnew(l,4);
        Qnew(l,N-1:N)=Qnew(l,N-2);
    end
    Qold=Qnew;
    for i=1:N
        Unew(1,i)=Qnew(2,i)/Qold(1,i);
    end
end
figure(1);
plot(Qnew(1,:));
title('Plot for Flow of Water and Height Variation using WENO');
xlabel('x in m'); ylabel('H in m');
pause(0.001);
figure(2);
plot(Unew(1,:));
title('Plot for Velocity of Water using WENO');
xlabel('x in m'); ylabel('U in m/s');
end

function [df]= Weno(Qold,delta_t)
N=101;
g=9.81;
for i=1:N
    cold(i)=sqrt(g.*Qold(1,i));
end
%Calculating eigen values
for i=1:N
        eigen(1,i)=(Qold(2,i)./Qold(1,i))+cold(i);
        eigen(2,i)=(Qold(2,i)./Qold(1,i))-cold(i);
end
beta=max(abs(eigen));
for i=1:N
    E(1,i)=Qold(2,i);
    E(2,i)=((Qold(2,i)^2/Qold(1,i))+(0.5*g*(Qold(1,i)^2)));
end
weight0=1/10; weight1=6/10; weight2=3/10;
for l=1:2
    for i=3:N-2
        b0(l,i)=13/12*(E(l,i-2)-2*E(l,i-1)+E(l,i))^2+1/4*(3*E(l,i-2)-4*E(l,i-1)+E(l,i))^2;
        b1(l,i)=13/12*(E(l,i-1)-2*E(l,i)+E(l,i+1))^2+1/4*(3*E(l,i-1)-E(l,i+1))^2;
        b2(l,i)=13/12*(E(l,i)-2*E(l,i+1)+E(l,i+2))^2+1/4*(E(l,i)-4*E(l,i+1)+E(l,i+2))^2;
    end
end
ep=10^-6;
for l=1:2
    for i=3:N-2
        wir0(l,i)=weight0/((ep+b0(l,i))^2+(ep+b1(l,i))^2+(ep+b2(l,i))^2);
        wir1(l,i)=weight1/((ep+b0(l,i))^2+(ep+b1(l,i))^2+(ep+b2(l,i))^2);
        wir2(l,i)=weight2/((ep+b0(l,i))^2+(ep+b1(l,i))^2+(ep+b2(l,i))^2);
    end
end
for l=1:2
    for i=3:N-2
        wzero(l,i)=wir0(l,i)/(wir0(l,i)+wir1(l,i)+wir2(l,i));
        wone(l,i)=wir1(l,i)/(wir0(l,i)+wir1(l,i)+wir2(l,i));
        wtwo(l,i)=wir2(l,i)/(wir0(l,i)+wir1(l,i)+wir2(l,i));
    end
end
for l=1:2
    for i=3:N-2
        s0(l,i)=1/3*E(l,i-2)-7/6*E(l,i-1)+11/6*E(l,i);
        s1(l,i)=-1/6*E(l,i-1)+5/6*E(l,i)+1/3*E(l,i+1);
        s2(l,i)=1/3*E(l,i)+5/6*E(l,i+1)-1/6*E(l,i+2);
    end
end
for l=1:2
    for i=3:N-2
        Epos(l,i)=wzero(l,i)*s0(l,i)+wone(l,i)*s1(l,i)+wtwo(l,i)*s2(l,i);
    end
end
for l=1:2
    for i=4:N-1
        c0(l,i)=13/12*(E(l,i-3)-2*E(l,i-2)+E(l,i-1))^2+1/4*(3*E(l,i-3)-4*E(l,i-2)+E(l,i-1))^2;
        c1(l,i)=13/12*(E(l,i-2)-2*E(l,i-1)+E(l,i))^2+1/4*(3*E(l,i-2)-E(l,i))^2;
        c2(l,i)=13/12*(E(l,i-1)-2*E(l,i)+E(l,i+1))^2+1/4*(E(l,i-1)-4*E(l,i)+E(l,i+1))^2;
    end
end
for l=1:2
    for i=4:N-1
        vir0(l,i)=weight0/((ep+c2(l,i))^2+(ep+c1(l,i))^2+(ep+c0(l,i))^2);
        vir1(l,i)=weight1/((ep+c2(l,i))^2+(ep+c1(l,i))^2+(ep+c0(l,i))^2);
        vir2(l,i)=weight2/((ep+c2(l,i))^2+(ep+c1(l,i))^2+(ep+c0(l,i))^2);
    end
end
for l=1:2
    for i=4:N-1
        vzero(l,i)=vir0(l,i)/(vir0(l,i)+vir1(l,i)+vir2(l,i));
        vone(l,i)=vir1(l,i)/(vir0(l,i)+vir1(l,i)+vir2(l,i));
        vtwo(l,i)=vir2(l,i)/(vir0(l,i)+vir1(l,i)+vir2(l,i));
    end
end
for l=1:2
    for i=4:N-1
        t0(l,i)=1/3*E(l,i-3)-7/6*E(l,i-2)+11/6*E(l,i-1);
        t1(l,i)=-1/6*E(l,i-2)+5/6*E(l,i-1)+1/3*E(l,i);
        t2(l,i)=1/3*E(l,i-1)+5/6*E(l,i)-1/6*E(l,i+1);
    end
end
for l=1:2
    for i=4:N-1
        Eneg(l,i)=vzero(l,i)*t0(l,i)+vone(l,i)*t1(l,i)+vtwo(l,i)*t2(l,i);
    end
end
for i=3:N-2
    Fplushalf(1,i)=0.5*(Epos(1,i)+Eneg(1,i+1))-0.5*max(beta(i),beta(i+1))*(Qold(1,i+1)-Qold(1,i));
    Fplushalf(2,i)=0.5*(Epos(2,i)+Eneg(2,i+1))-0.5*max(beta(i),beta(i+1))*(Qold(2,i+1)-Qold(2,i));
   % Fplushalf(3,i)=0.5*(Epos(3,i)+Eneg(3,i+1))-0.5*max(beta(i),beta(i+1))*(Qold(3,i+1)-Qold(3,i));
end
for i=4:N-1
    Fminhalf(1,i)=0.5*(Eneg(1,i)+Epos(1,i-1))-0.5*max(beta(i),beta(i-1))*(Qold(1,i)-Qold(1,i-1));
    Fminhalf(2,i)=0.5*(Eneg(2,i)+Epos(2,i-1))-0.5*max(beta(i),beta(i-1))*(Qold(2,i)-Qold(2,i-1));
    %Fminhalf(3,i)=0.5*(Eneg(3,i)+Epos(3,i-1))-0.5*max(beta(i),beta(i-1))*(Qold(3,i)-Qold(3,i-1));
end
for k=1:2
    for i=4:N-2
        df(k,i)=(Fplushalf(k,i)-Fminhalf(k,i));
    end
end
end