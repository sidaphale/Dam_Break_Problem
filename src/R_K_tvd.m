%Runge-Kutta with TVD Scheme
close all; clear all; clc;
N=101; 
L=100; 
dx=L/(N-1); 
g=9.81; 
dt=0.001; 
t=4; 
c=dt/dx; 
% intial values of h and v
for i=1:N
    if i*dx<50
        h(i)=10;
        v(i)=0;
    else
        h(i)=1;
        v(i)=0;
    end
end
% intializing variables for H_old and H_new to zero H_old is for h and
% U_old is for v
for i=1:N
    H_old(i)=0;
    U_old(i)=0;
    H_new(i)=0;
    U_new(i)=0;
    alphah1(i)=0;
    alphah2(i)=0;
    alphav1(i)=0;
    alphav2(i)=0;
    sigmah1(i)=0;
    sigmah2(i)=0;
    sigmav1(i)=0;
    sigmav2(i)=0;
    rh1(i)=0;
    rh2(i)=0;
    rv1(i)=0;
    rv2(i)=0;
    Gh1(i)=0; Gh2(i)=0; Gv1(i)=0; Gv2(i)=0;
    phih1(i)=0; phih2(i)=0; phiv1(i)=0; phiv2(i)=0;
end
% intializing the U_old values to their corresponding h and v
for i=1:N
    H_old(i)=h(i);
    U_old(i)=v(i)*h(i);
end
for i=1:N
    H1(i)=0; H2(i)=0; H3(i)=0; H4(i)=0; %Defining R-K terms for height
    U1(i)=0; U2(i)=0; U3(i)=0; U4(i)=0; %Defining R-K terms for velocity
end
tinitial=0;
while tinitial<=t
    tinitial=tinitial+dt;
    %Defining Roe-Sweby Flux Limiter
    for i=1:N-1
        if ((H_old(i+1)-H_old(i))>=-0.1) || ((H_old(i+1)-H_old(i))<=0.1)
            alphah1(i)=H_old(i);
        else
            alphah1(i)=(U_old(i+1)-U_old(i))/(H_old(i+1)-H_old(i));
        end
        if ((U_old(i+1)-U_old(i))>=-0.1) || ((U_old(i+1)-U_old(i))<=0.1)
            alphav1(i)=U_old(i);
        else
            alphav1(i)=((((U_old(i+1).^2)/H_old(i+1))+0.5*g*H_old(i+1).^2)-((U_old(i).^2)./H_old(i)+0.5*g*H_old(i).^2)./(U_old(i+1)-U_old(i)));
        end
    end
    for i=2:N-1
        if ((H_old(i+1)-H_old(i))>=-0.1) || ((H_old(i+1)-H_old(i))<=0.1)
            alphah2(i)=H_old(i);
        else
            alphah2(i)=(U_old(i)-U_old(i-1))/(H_old(i)-H_old(i-1));
        end
        if ((U_old(i+1)-U_old(i))>=-0.1) || ((U_old(i+1)-U_old(i))<=0.1)
            alphav2(i)=U_old(i);
        else
            alphav2(i)=((((U_old(i).^2)/H_old(i)+0.5*g*H_old(i).^2)-((U_old(i-1).^2)./H_old(i-1)+0.5*g*H_old(i-1).^2))./(U_old(i)-U_old(i-1)));
        end
    end
    %Finding values of sigma
    for i=1:N
        if alphah1(i)==0
            sigmah1(i)=0;
        else
            sigmah1(i)=alphah1(i)/abs(alphah1(i));
        end
    end
    for i=1:N
        if alphav1(i)==0
            sigmav1(i)=0;
        else
            sigmav1(i)=alphav1(i)/abs(alphav1(i));
        end
    end
    for i=2:N
        if alphah2(i)==0
            sigmah2(i)=0;
        else
            sigmah2(i)=alphah2(i)/abs(alphah2(i));
        end
    end
    for i=2:N
        if alphav2(i)==0
            sigmav2(i)=0;
        else
            sigmav2(i)=alphav2(i)/abs(alphav2(i));
        end
    end
    %Calculating the values of r
    for i=2:N-2
        if (H_old(i+1)-H_old(i)==0)
            rh1(i)=0;
        else
            rh1(i)=(H_old(i+1+sigmah1(i))-H_old(i+sigmah1(i)))/(H_old(i+1)-H_old(i));
        end
    end
    for i=2:N-2
        if (U_old(i+1)-U_old(i)==0)
            rv1(i)=0;
        else
            rv1(i)=(U_old(1+i+sigmav1(i))-U_old(i+sigmav1(i)))/(U_old(i+1)-U_old(i));
        end
    end
    for i=3:N-2
        if (H_old(i)-H_old(i-1)==0)
            rh2(i)=0;
        else
            rh2(i)=(H_old(i+sigmah2(i))-H_old(i-1+sigmah2(i)))/(H_old(i)-H_old(i-1));
        end
    end
    for i=3:N-2
        if (U_old(i)-U_old(i-1)==0)
            rv2(i)=0;
        else
            rv2(i)=(U_old(i+sigmav2(i))-U_old(i-1+sigmav2(i)))/(U_old(i)-U_old(i-1));
        end
    end
    %Calculate values of G
    for i=3:N-2
        Gh1(i)=(rh1(i)+abs(rh1(i)))/(1+rh1(i));
    end
    for i=3:N-2
        Gh2(i)=(rh2(i)+abs(rh2(i)))/(1+rh2(i));
    end
    for i=3:N-2
        Gv1(i)=(rv1(i)+abs(rv1(i)))/(1+rv1(i));
    end
    for i=3:N-2
        Gv2(i)=(rv2(i)+abs(rv2(i)))/(1+rv2(i));
    end   
    %Calculate the values of Flux Limiter phi(i+1/2) and phi(i-1/2)
    for i=3:N-2
        phih1(i)=((Gh1(i)/2)*(abs(alphah1(i))+c*(alphah1(i).^2))-alphah1(i))*(H_old(i+1)-H_old(i));
    end
    for i=3:N-2
        phih2(i)=((Gh2(i)/2)*(abs(alphah2(i))+c*(alphah2(i).^2))-alphah2(i))*(H_old(i)-H_old(i-1));
    end
    for i=3:N-2
        phiv1(i)=((Gv1(i)/2)*(abs(alphav1(i))+c*(alphav1(i).^2))-alphav1(i))*(U_old(i+1)-U_old(i));
    end
    for i=3:N-2
        phiv2(i)=((Gv2(i)/2)*(abs(alphav2(i))+c*(alphav2(i).^2))-alphav2(i))*(U_old(i)-U_old(i-1));
    end
    %Implementing Runge-Kutta Method
    for i=1:N
        H1(i)=H_old(i); U1(i)=U_old(i);
    end
    for i=2:N-1
        H2(i)=H_old(i)-((0.25*c*0.5.*(H1(i+1).*U1(i+1)-H1(i-1).*U1(i-1)))/(H1(i)));
        U2(i)=U_old(i)-(0.25*c*0.5.*(((U1(i+1).^2/H1(i+1))+0.5*g*(H1(i+1).^2))-((U1(i-1).^2./H1(i-1))+0.5*g*H1(i-1).^2)));
    end
    H2(1)=H_old(1); U2(1)=U_old(1);
    H2(N)=H_old(N-1); U2(N)=U_old(N-1);
    for i=2:N-1
        H3(i)=H_old(i)-(((1/3)*c*0.5.*(H2(i+1).*U2(i+1)-H2(i-1).*U2(i-1)))./(H2(i)));
        U3(i)=U_old(i)-((1/3)*c*0.5.*(((U2(i+1).^2./H2(i+1))+0.5*g*H2(i+1).^2)-((U2(i-1).^2./H2(i-1))+0.5*g*H2(i-1).^2)));
    end
    H3(1)=H_old(1); U3(1)=U_old(1);
    H3(N)=H_old(N-1); U3(N)=U_old(N-1);
    for i=2:N-1
        H4(i)=H_old(i)-((0.5*c*0.5.*(H3(i+1).*U3(i+1)-H3(i-1).*U3(i-1)))./(H3(i)));
        U4(i)=U_old(i)-(0.5*c*0.5.*(((U3(i+1).^2./H3(i+1))+0.5*g*H3(i+1).^2)-((U3(i-1).^2./H3(i-1))+0.5*g*H3(i-1).^2)));
    end
    H4(1)=H_old(1); U4(1)=U_old(1);
    H4(N)=H_old(N-1); U4(N)=U_old(N-1);
    for i=2:N-1
        H_new(i)=H_old(i)-0.25*(c*0.5.*(H4(i+1).*U4(i+1)-H4(i-1).*U4(i-1))./(H4(i)));
        U_new(i)=U_old(i)-0.25*(c*0.5.*(((U4(i+1).^2./H4(i+1))+0.5*g*H4(i+1).^2)-((U4(i-1).^2./H4(i-1))+0.5*g*H4(i-1).^2)));
    end
    H_new(1)=H_old(1); U_new(1)=U_old(1);
    H_new(N)=H_old(N-1); U_new(N)=U_old(N-1);
    %Adding the flux limiter
    for i=1:N
        H_new(i)=H_new(i)-0.5*c.*(phih1(i)-phih2(i));
        U_new(i)=U_new(i)-0.5*c.*(phiv1(i)-phiv2(i));
    end
    for i=1:N
        H_old(i)=H_new(i);
        U_old(i)=U_new(i);
    end
    Unewfinal=U_new./H_new;
figure(1);
plot(H_new);
title('Plot for water flow by Runge-Kutta with TVD');
xlabel('x in m'); ylabel('H in m');
pause(0.001);
end
figure(1);
plot(H_new);
title('Plot for water flow by Runge-Kutta with TVD');
xlabel('x in m'); ylabel('H in m');
figure(2);
plot(Unewfinal);
title('Plot for velocity by Runge-Kutta with TVD');
xlabel('x in m'); ylabel('Velocity U in m/s'); 