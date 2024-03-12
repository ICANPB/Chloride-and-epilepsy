G_NaL=0.04;G_KL=0.1;G_ClL=0.15;gamma=0.03;
belta=1;g_Na=30;g_K=20;tao=1000;
G_gaba=0;rou=0.2;HCO_o=26;HCO_i=12;

c_K_t=100;c_Na_t=155;c_Cl_t=151;
l_t=50000000;dt=0.1;stepc=1;t_tsee=stepc:stepc:l_t;l_tc=length(t_tsee);
Kbath=20;ipc_all=[0.25,0.025,0.0025];lll=length(ipc_all);

X_shif_f2=zeros(l_tc,7,1);

rate=5;tao_d=10;tao_r=1;t_all=0:dt:3000;a=0;
feq_gaba=10;di=floor((1000/feq_gaba)/dt);inhib_sp_po=zeros(1,floor(l_t/di));
f_gaba=rate*(exp(-t_all/tao_d)-exp(-t_all/tao_r))/(tao_d-tao_r);

for i2=1:1
    
    aaa=rem(i2,3);
    if aaa==0
        aaa=3;
    end
    bbb=(i2-aaa)/3+1;
    
    X_RTM=X_shif_k3(1:7);
    
    G_RTM=[G_gaba,rou,HCO_o,HCO_i];
    
    Y_RTM=[g_K,G_KL,g_Na,G_NaL,G_ClL,c_Na_t,c_Cl_t,gamma,belta,tao,dA];
    
    
    for i1=1:l_t
        
        tt1=0.1*l_t;tt2=l_t;
        if i1>tt1
            if rem(i1,di)==0
                a=a+1;
                inhib_sp_po(a)=i1-1;
            end
            if a>0
                aa=i1-inhib_sp_po(1:a);
                aa(aa>3000)=[];
                G_gaba=sum(f_gaba(aa));
            end
        end
        if i1>tt2
            G_gaba=0;
        end
        G_RTM(1)=G_gaba;
        k1=RTM_gaba(X_RTM,Y_RTM,G_RTM);k2=RTM_gaba(X_RTM+(dt/2)*k1,Y_RTM,G_RTM);k3=RTM_gaba(X_RTM+(dt/2)*k2,Y_RTM,G_RTM);k4=RTM_gaba(X_RTM+dt*k3,Y_RTM,G_RTM);
        X_RTM(1:7)=X_RTM(1:7)+(dt/6)*(k1(1:7)+2*k2(1:7)+2*k3(1:7)+k4(1:7));%X_RTM(9)=k1(9);
        
        if rem(i1,stepc)==0
            ii1=i1/stepc;
            X_shif_f2(ii1,:,i2)=X_RTM;
        end
    end
end