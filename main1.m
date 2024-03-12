Glu1=0;Glu=0;G_NaL=0.04;G_KL=0.1;G_ClL=0.1;gamma=0.03;
belta=7;g_Na=30;g_K=20;tao=1000;
G_gaba=0;rou=0;HCO_ob=24;HCO_ib=14.7;HCO_o=24;HCO_i=14.7;
ipc_hi=0;ipc_ho=0;U_KCC2=0.3;U_NKCC1=0.1;

c_K_t=(100-3)+3/belta;c_Na_t=275;c_Cl_t=187;
l_t=4000000;dt=0.1;stepc=1;t_tsee=stepc:stepc:l_t;l_tc=length(t_tsee);
Kbath=20;ipc=0;Kbath_all=linspace(3,15,30);
G_RTM=[G_gaba,rou,HCO_ob,HCO_ib,ipc_ho,ipc_hi];

X_reco=zeros(l_tc,10,1);Gex=zeros(l_tc,2);

rate=5;tao_d=10;tao_r=1;t_all=0:dt:3000;a=0;
feq_gaba=80;di=floor((1000/feq_gaba)/dt);inhib_sp_po=zeros(1,floor(l_t/di));
f_gaba=rate*(exp(-t_all/tao_d)-exp(-t_all/tao_r))/(tao_d-tao_r);
kkk_all=[28,28];kk=16.47:0.01:16.57;

for i2=1:4

    cli=[3.89,3.90,6.02,6.03];
    X_RTM=X_k3_b7_cv;X_RTM(7)=X_RTM(7)+cli(i2);

    K_diff=[3,0];
    V=X_RTM(1);c_K_o=X_RTM(5);c_Cl_i=X_RTM(7);c_Na_i=X_RTM(6);
    dA=84.4695;
    Y_RTM=[g_K,G_KL,g_Na,G_NaL,G_ClL,c_Na_t,c_Cl_t,gamma,belta,tao,dA,Glu,U_KCC2,U_NKCC1];

    rate=5;tao_d=10;tao_r=1;t_all=0:dt:3000;a=0;
    feq_gaba=80;di=floor((1000/feq_gaba)/dt);inhib_sp_po=zeros(1,floor(l_t/di));
    f_gaba=rate*(exp(-t_all/tao_d)-exp(-t_all/tao_r))/(tao_d-tao_r);
    
    for i1=1:l_t
        
        tt1=l_t;tt2=0.8*l_t;
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
        else
            G_gaba=0;
        end

        G_RTM(1)=G_gaba;
        
        k1=RTM_Kdiffbi_cv(X_RTM,Y_RTM,G_RTM,K_diff);k2=RTM_Kdiffbi_cv(X_RTM+(dt/2)*k1,Y_RTM,G_RTM,K_diff);
        k3=RTM_Kdiffbi_cv(X_RTM+(dt/2)*k2,Y_RTM,G_RTM,K_diff);k4=RTM_Kdiffbi_cv(X_RTM+dt*k3,Y_RTM,G_RTM,K_diff);
        X_RTM=X_RTM+(dt/6)*(k1+2*k2+2*k3+k4);
        
        for i3=2:4
            if X_RTM(i3)>1
                X_RTM(i3)=1;
            elseif X_RTM(i3)<0
                X_RTM(i3)=0;
            end
        end

        if rem(i1,stepc)==0
            ii1=i1/stepc;
            X_reco(ii1,:,i2)=X_RTM;
            Gex(ii1,1,i2)=Glu;Gex(ii1,2,i2)=G_gaba;
        end
    end

end