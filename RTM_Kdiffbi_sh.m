function [XX] = RTM_Kdiffbi_sh(X,Y,G,K)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
V=X(1);m=X(2);h=X(3);n=X(4);c_K_o=X(5);c_Na_i=X(6);c_Cl_i=X(7);c_K_i=X(8);HCO_o=X(9);HCO_i=X(10);
g_K=Y(1);G_KL=Y(2);g_Na=Y(3);G_NaL=Y(4);G_ClL=Y(5);c_Na_t=Y(6);c_Cl_t=Y(7);gamma=Y(8);belta=Y(9);tao=Y(10);dA=Y(11);Glu=Y(12);
%U_KCC2=Y(13);U_NKCC1=Y(14);
Kbath=K(1);ipc=K(2);
c_K_i=2*((gamma/1000)*V-((-c_Na_t+c_Cl_t)/(2*belta)+c_Na_i-c_K_o/(2*belta)-c_Cl_i)-dA);
c_K_i=(gamma/1000)*V+c_Cl_i-c_Na_i+dA;
G_gaba=G(1);rou=G(2);HCO_ob=G(3);HCO_ib=G(4);ipc_Ho=G(5);ipc_Hi=G(6);
c_Na_o=c_Na_t-belta*c_Na_i;c_Cl_o=c_Cl_t-belta*c_Cl_i;
I_pump=(0.8/(1+exp((25-c_Na_i)/3)))*(1/(1+exp(3.5-c_K_o)));
U_KCC2=0.3;
I_KCC2=U_KCC2*log((c_K_i*c_Cl_i)/(c_K_o*c_Cl_o));
if I_KCC2<0
    I_KCC2=0;
end
U_NKCC1=0.1;
I_NKCC1=(U_NKCC1/(1+exp(16-c_K_o)))*(log((c_K_i*c_Cl_i)/(c_K_o*c_Cl_o))+log((c_Na_i*c_Cl_i)/(c_Na_o*c_Cl_o)));
%I_NKCC1=(U_NKCC1/(1))*(log((c_K_i*c_Cl_i)/(c_K_o*c_Cl_o))+log((c_Na_i*c_Cl_i)/(c_Na_o*c_Cl_o)));
if I_NKCC1>0
    I_NKCC1=0;
end
%I_KCC2=0;I_NKCC1=0;
%V=((-c_Na_t+c_Cl_t)/2+c_Na_i+(c_K_i-c_K_o)/2-c_Cl_i)*(1000/gamma)+dA;
E_Na=26.64*log(c_Na_o/c_Na_i);E_K=26.64*log(c_K_o/c_K_i);E_Cl=26.64*log(c_Cl_i/c_Cl_o);
E_HCO=26.64*log(HCO_i/HCO_o);
XX(1)=-g_Na*(m^3)*h*(V-E_Na)-g_K*(n^4)*(V-E_K)-(G_NaL+0.4*Glu)*(V-E_Na)-(G_KL+0.2*Glu)*(V-E_K)-(G_ClL+G_gaba)*(V-E_Cl)-rou*G_gaba*(V-E_HCO)-I_pump/gamma;
XX(2)=((0.32*(V+54))/(1-exp(-(V+54)/4)))*(1-m)-((0.28*(V+27))/(exp((V+27)/5)-1))*m;
XX(3)=(0.128*exp(-(V+50)/18))*(1-h)-(4/(1+exp(-(V+27)/5)))*h;
XX(4)=((0.032*(V+52))/(1-exp(-(V+52)/5)))*(1-n)-(0.5*exp(-(V+57)/40))*n;
XX(5)=(((g_K*(n^4)*(V-E_K)+(G_KL+0.2*Glu)*(V-E_K))*gamma-2*I_pump+I_NKCC1+I_KCC2)*belta-ipc*(c_K_o-Kbath))/tao;
XX(6)=((-g_Na*(m^3)*h*(V-E_Na)-(G_NaL+0.4*Glu)*(V-E_Na))*gamma-3*I_pump-I_NKCC1)/tao;
XX(7)=(((G_ClL+G_gaba)*(V-E_Cl))*gamma-2*I_NKCC1-I_KCC2)/tao;
XX(8)=((-(g_K*(n^4)*(V-E_K)+(G_KL+0.2*Glu)*(V-E_K))*gamma+2*I_pump-I_NKCC1-I_KCC2))/tao;
XX(9)=(-(rou*G_gaba*(V-E_HCO))*gamma*belta-ipc_Ho*(HCO_o-HCO_ob))/tao;
XX(10)=((rou*G_gaba*(V-E_HCO))*gamma-ipc_Hi*(HCO_i-HCO_ib))/tao;
end