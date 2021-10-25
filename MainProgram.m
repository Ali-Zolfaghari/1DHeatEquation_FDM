
%***************************************************************************************************
%*   Solve 1D-heat equation by presented code.
%*   I take no responsibilities for any errors in the code or damage thereby.
%*   Please notify me at zolfaghari1992iut@gmail.com if the code is used in any type of application.
%***************************************************************************************************
%*   Developer   : Ali Zolfaghari Sichani (10-03-2017)
%***************************************************************************************************
%*   References  : 
%*   Computational Fluid Mechanics and Heat Transfer.
%*   by John C. Tannehill (Author), Dale Anderson (Author), Richard H. Pletcher (Author).
%***************************************************************************************************
%*   Heat Equation in one-dimensional domain by Dirichlet boundary condition. (solving by finite differences methods)   :   
%*   Ut = a . Uxx
%*   Inputs      :
%*   Methods ID
%*               SimpleExplicit = 1
%*               SimpleImplicit = 2
%*                CrankNicolson = 3
%*                 CombineTypeA = 4
%*                 CombineTypeB = 5
%*                DufortFrankel = 6
%*   L          (domain length                                            )
%*   M          (number of division of domain                             )
%*   dt         (time step                                                )
%*   FinalTime  (stop criteria time                                       )
%*   Tl         (temperature of left boundary                             )
%*   Tr         (temperature of right boundary                            )
%*   T0         (initial temperature                                      )
%*   MAXERROR   (max. allowable error                                     )
%*   alpha      (underrelaxtion factor                                    )
%*   Solver     (1: MATLAB solver - 2: Successive over-relaxation solver  )
%*   Outputs     :
%*   plot numerical and exact solution of temperature distribution
%***************************************************************************************************


clear,clc,close all
format compact
format long


%% input
L = 1.0;
M = 20;
dt = 0.001;
FinalTime = 0.5;
alpha = 1.0;
Tl = 0.0;
Tr = 0.0;
T0 = 100.0;
MAXERROR = 0.001;
Ws = 0.85;
Solver = 1;



%% initial
S = {'SimpleExplicit','SimpleImplicit','CrankNicolson','CombineTypeA','CombineTypeB','DufortFrankel'};
dx = L/M;
x = 0.0:dx:L;
Teta_A = 0.5-((dx*dx)/(12.0*alpha*dt));
Teta_B = 0.5+((dx*dx)/(12.0*alpha*dt));

%% exact
T_ext = zeros(1,M+1);
for i = 1:M+1
    for n = 1:1000
        K = n*pi;
        An = -(2*T0*(cos(K) - 1))/K;
        T_ext(i) = T_ext(i)+An*exp(-alpha*K*K*FinalTime)*sin(K*x(i));
    end
end

%% numerical
for IMethod = 1:6
    
    T_old = T0*ones(1,M+1);
    T_new = T0*ones(1,M+1);
    time = 0.0;
    
    if(IMethod == 1)
        
        while(time < FinalTime)
            [ A,B ] = SimpleExplicit( T_old,x,dx,dt,M,Tl,Tr,alpha );
            if(Solver == 1)
                U = A\B;
                T_new = U';
            else
                [ T_new ] = SOR( A,B,Ws,MAXERROR );
            end
            T_old = T_new;
            time = time+dt;
        end
        Tf = T_new;
        
    elseif(IMethod == 2)
        
        while(time < FinalTime)
            [ A,B ] = SimpleImplicit( T_old,x,dx,dt,M,Tl,Tr,alpha );
            if(Solver == 1)
                U = A\B;
                T_new = U';
            else
                [ T_new ] = SOR( A,B,Ws,MAXERROR );
            end
            T_old = T_new;
            time = time+dt;
        end
        Tf = T_new;
        
    elseif(IMethod == 3)
        
        while(time < FinalTime)
            [ A,B ] = CrankNicolson( T_old,x,dx,dt,M,Tl,Tr,alpha );
            if(Solver == 1)
                U = A\B;
                T_new = U';
            else
                [ T_new ] = SOR( A,B,Ws,MAXERROR );
            end
            T_old = T_new;
            time = time+dt;
        end
        Tf = T_new;
        
    elseif(IMethod == 4)
        
        while(time < FinalTime)
            [ A,B ] = CombineTypeA( T_old,x,dx,dt,M,Tl,Tr,alpha,Teta_A );
            if(Solver == 1)
                U = A\B;
                T_new = U';
            else
                [ T_new ] = SOR( A,B,Ws,MAXERROR );
            end
            T_old = T_new;
            time = time+dt;
        end
        Tf = T_new;
        
    elseif(IMethod == 5)
        
        [ A,B ] = SimpleExplicit( T_old,x,dx,dt,M,Tl,Tr,alpha );
        if(Solver == 1)
            U = A\B;
            T_new = U';
        else
            [ T_new ] = SOR( A,B,Ws,MAXERROR );
        end
        T2 = T_new;
        T1 = T_old;
        time = dt;
        while(time < FinalTime)
            [ A,B ] = CombineTypeB( T1,T2,x,dx,dt,M,Tl,Tr,alpha,Teta_B );
            if(Solver == 1)
                U = A\B;
                T_new = U';
            else
                [ T_new ] = SOR( A,B,Ws,MAXERROR );
            end
            T1 = T2;
            T2 = T_new;
            time = time+dt;
        end
        Tf = T_new;
        
    elseif(IMethod == 6)
        
        [ A,B ] = SimpleExplicit( T_old,x,dx,dt,M,Tl,Tr,alpha );
        if(Solver == 1)
            U = A\B;
            T_new = U';
        else
            [ T_new ] = SOR( A,B,Ws,MAXERROR );
        end
        T2 = T_new;
        T1 = T_old;
        time = dt;
        while(time < FinalTime)
            [ A,B ] = DufortFrankel( T1,T2,x,dx,dt,M,Tl,Tr,alpha );
            if(Solver == 1)
                U = A\B;
                T_new = U';
            else
                [ T_new ] = SOR( A,B,Ws,MAXERROR );
            end
            T1 = T2;
            T2 = T_new;
            time = time+dt;
        end
        Tf = T_new;
        
    end
    
    T{IMethod} = Tf;
    Error = (norm(Tf-T_ext))/(M+1);
    fig = figure(IMethod);hold on;grid on;
    plot(x,T_ext,'b-','linewidth',2.2);
    plot(x,Tf,'r-*','linewidth',1.1);
    xlabel('X','FontWeight','bold','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
    ylabel('Temperature','FontWeight','bold','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
    title(['Time = ',num2str(FinalTime),', Error = ',num2str(Error)],'FontWeight','bold','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
    legend({'Exact',S{IMethod}});
    saveas(fig,[num2str(IMethod),'.jpg']);
    
end

cl = jet(7);
fig = figure(7);hold on;grid on;
plot(x,T_ext,'color',cl(1,:),'linewidth',2.5);
for IMethod = 1:6
    plot(x,T{IMethod},'color',cl(IMethod+1,:),'linewidth',1.25);
end
xlabel('X','FontWeight','bold','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
ylabel('Temperature','FontWeight','bold','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
title(['Time = ',num2str(FinalTime)],'FontWeight','bold','FontSize',12,'FontName','Times New Roman','FontAngle','italic');
legend({'Exact','SimpleExplicit','SimpleImplicit','CrankNicolson','CombineTypeA','CombineTypeB','DufortFrankel'});
saveas(fig,'7.jpg');

