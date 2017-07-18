% Gina M. Eberhart
% GE BEMT Thrust Prediction
clear all;
clc;
close all;
%% Initial Conditions and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import Tabulated Airfoil Data
clarky = Airfoilimport('Clarky 100k.txt.txt');
% Load Twist and Taper Provile for Propeller
load TwistandTaperData10x4.mat;
%Pitch (rads)
th=deg2rad(theta);
% Radial Chord (m)
chord=cinterp;
% Rotor Radius (m)
R=.123;
% Area of Rotor Disk
A=pi*R^2;
% Number of Propellers
NP=4;
% Acceleration Due to Gravity (m/s^2)
g=9.81;
% Mass of Payload(kg)
MP=0;
% Payload Weight (N)
WP=MP*g;
% Totall Mass of Quadcopter(kg)
MT=1.435;
% Total Weight of Quadcopter (N)
WT=MT*g+WP;
% Thrust Required to Hover From Single Propeller (N)
ThrustSP=WT/NP;
% Density of Air (kg/m^3)
rho=1.225;
% Number of Elements to Analyze
N=20;
% Length of Element (m)
dr=R/N;
% Radial Distance of Elements
r1=[dr:dr:R];
% 10x4.7 thrust curve
c1=0.2685;
c2=-1.312;
c3=0.9949;
c4=0.001212;
% Tip Speed Ratio
TSR=41.35;
%Height Ratio
HR=[0.5:0.1:6];
%Rotor RPM
RPM=[4000:250:10000];
%% Calculation of Power Requirements with Consideration for GE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for l=1:length(HR)
    for k=1:length(RPM)
        omega(k)=(RPM(k)/60)*2*pi;
        for i=1
            % Free stream velocity guess
            V(i)=(omega(k)*R)/TSR(i);
            % Initialize thrust and torque
            T=0;
            Q=0;
            for j=1:length(r1)
                % Going Along the Radius of the blade rad=radius
                r=r1(j);
                % Initial Guess for Axial Induction
                a=0.1;
                % Initial Guess for Tangential Induction
                b=0.01;
                % Setting a value for while loop conditions
                TerminateSEQ=0;
                count=1;
                while (TerminateSEQ==0),
                    % Solving for V0 and V2 using free stream velocity and induction factor
                    V0=V(i)*(1+a);
                    V2=omega(k)*r*(1-b);
                    % Solving for phi using velocity components
                    phi=atan2(V0,V2);
                    % Solving for angle of attack
                    alpha=th(j)-phi;
                    %Finding CL and CD based on angle of attack from Tabulated Data
                    if rad2deg(alpha)>=min(clarky.alpha) && rad2deg(alpha)<max(clarky.alpha)
                        cl=interp1(deg2rad(clarky.alpha),clarky.CL,alpha);
                        cd=interp1(deg2rad(clarky.alpha),clarky.CD,alpha);
                    else
                        cl=0;
                        cd=1;
                    end
                    % Finding local velocity
                    Vlocal=sqrt(V0^2+V2^2);
                    %Coefficients of thrust and torque
                    ctio=c1*exp(c2*HR(l))+c3*exp(c4*HR(l));
                    ct=(cl*cos(phi)-cd*sin(phi))*ctio;
                    cq=(cd*cos(phi)+cl*sin(phi));
                    
                    % Finding Differential Thrust BET
                    DtBET=0.5*rho*Vlocal^2*2*chord(j)*ct;
                    
                    % Finding Differential Torque BET
                    DqBET=0.5*rho*Vlocal^2*2*chord(j)*r*cq;
                    
                    % Finding Differential Thrust MT
                    DtMT=(4*pi*r*rho*V(i)^2*(1+a));
                    % Finding Differential Torque MT
                    DqMT=(4*pi*r^3*rho*V(i)*(1+a)*omega(k));
                    % BEMT thrust comparison
                    Tratio=DtBET/DtMT;
                    Qratio=DqBET/DqMT;
                    % Cal
                    anew=0.5*(a+Tratio);
                    bnew=0.5*(b+Qratio);
                    if (abs(anew-a)<1.0e-5),
                        if (abs(bnew-b)<1.0e-5),
                            TerminateSEQ=1;
                        end;
                    end;
                    % Reset induction factors
                    a=anew;
                    b=bnew;
                    count=count+1;
                    if (count>500),
                        TerminateSEQ=1;
                    end;
                end;
                phieach(j)=rad2deg(phi);
                cleach(j)=cl;
                alphadegs(j)=rad2deg(alpha);
                
                T=T+DtBET*dr;
                Q=Q+DqBET*dr;
            end;
            thrustis(k,i)=T;
            torqueis(k,i)=Q;
        end
    end    
    % Rotational Speed to produce 1/4 of T Required to Hover (for quadcopter)
    Omegareq(l)=interp1(thrustis,omega,ThrustSP);
    % Torque OGE for Single Prop
    TorqueOGE(l)=interp1(omega,torqueis,Omegareq(l));
    % Coefficeint of Torque for Single Prop
    CqOGE(l)=TorqueOGE(l)/(rho*A*Omegareq(l)^2*R^3);
    % Equating Coefficient of Power and Torque
    CpOGE(l)=CqOGE(l);
    % Power Required OGE for Single Prop
    P=rho*A*(Omegareq(l)*R)^3*CpOGE(l);
    %Motor Efficiency (70%)
    MEff=0.70;
    %Total Power Required
    Ptot=NP*P;
    Ptoteff(l)=Ptot+Ptot*MEff;    
    % Coefficient of Thrust
    Ct(l)=ThrustSP/(rho*(Omegareq(l)*R)^2*A);    
    % RPM Required
    RPMreq(l)=(Omegareq(l)*60)/(2*pi);
end

%% Figures
figure
plot(HR,Ft,'o')
xlabel('HR (z/R)')
ylabel('Flight Time (Mins)')
title('Flight Time vs. HR')
grid on

figure
plot(HR,Ptoteff,'o')
xlabel('HR (z/R)')
ylabel('Power IGE (Watts)')
title('Power vs. HR')
grid on

figure
plot(HR,RPMreq,'o')
xlabel('HR (z/R)')
ylabel('RPM(-)')
title('RPM vs. HR')
grid on

