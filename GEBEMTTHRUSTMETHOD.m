% Gina M. Eberhart
% GE BEMT Thrust Prediction
clear all;
clc;
close all;
%% Initial Conditions and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Airfoil Data Import (cl and cd vs. alpha)
clarky = Airfoilimport('Naca44121ookren.txt');
%Load measured 11x7 twist and taper profiles
load TwistandTaperData11x7.mat;
%Convert pitch from degrees to radians
th=deg2rad(theta);
%Chord Length (m)
chord=cinterp;
% Rotor Radius (m)
R=.139;
%Rotor Rotational Speed (RPM)
RPM=5000;
% Density of Air (standard)
rho=1.225;
% Converting Rotational Speed (RPM) to (rads/s)
omega=(RPM/60)*2*pi;
% Number of Elements
N=20;
%Size of Element (m)
dr=R/N;
% Span-wise Element Distances (m)
rs=[dr:dr:R];
% Tip Speed Ratio Range
TSR=[6:0.005:9];
% Height Ratio
HR=[0:.1:6];
% Altitude (m)
z=HR*R;
% Rotor Disk Area (m^2)
A=pi*R^2;
%Constants of Curve Fit
c1=0.0974;
c2=-1.935;
c3=1.005;
c4=-0.00049;

%% Calculation of Radial Thrust and Torque With Consideration for GE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1: length(TSR)
    % Free stream velocity guess
    V(i)=(omega*R)/TSR(i);
    for k=1:length(HR)
        % Initialize thrust and torque
        T=0;
        Q=0;
        for j=1:length(rs)
            %Going Along the Radius of Blade
            r=rs(j);
            % Initial Guess for Axial Induction
            a=0.1;
            %Initial Guess for Tangential Induction
            b=0.01;
            % Set termination condition initial value and count
            TerminateSEQ=0;
            count=1;
            while (TerminateSEQ==0),
                % Solving for V0 and V2 using free stream velocity and induction factor
                V0=V(i)*(1+a);
                V2=omega*r*(1-b);
                % Solving for phi using velocity components
                phi=atan2(V0,V2);
                % Solving for angle of attack
                alpha=th(j)-phi;
                % Finding CL and CD based on angle of attack from tabulated data
                if rad2deg(alpha)>=-9.5 && rad2deg(alpha)<17
                    cl=interp1(deg2rad(clarky.alpha),clarky.CL,alpha);
                    cd=interp1(deg2rad(clarky.alpha),clarky.CD,alpha);
                else
                    cl=0;
                    cd=1;
                end
                % Local velocity
                Vlocal=sqrt(V0^2+V2^2);
                % Modification of coefficient of thrust
                ctio(k,j)=c1*exp(c2*HR(k))+c3*exp(c4*HR(k));
                % Coefficient of thrust
                ct(k,j)=(cl*cos(phi)-cd*sin(phi));
                % Negate modification for non-thrust producing elements
                if ct(k,j)<0
                    ct(k,j)=ct(k,j);
                elseif ct(k,j)>0
                    ct(k,j)=(cl*cos(phi)-cd*sin(phi))*ctio(k,j);
                end
                % Coefficient of torque
                cq=(cd*cos(phi)+cl*sin(phi));
                % Finding differential BET thrust
                DtBET=0.5*rho*Vlocal^2*2*chord(j)*ct(k,j)*dr;
                % Finding differential BET torque
                DqBET=0.5*rho*Vlocal^2*2*chord(j)*r*cq*dr;
                % Finding differential MT thrust
                DtMT=4*pi*r*rho*V(i)^2*(1+a);T
                % Finding differential MT torque
                DqMT=4*pi*r^3*rho*V(i)*(1+a)*omega;
                % Comparison of BET to MT thrust
                Tratio=DtBET/DtMT;
                % Comparison of BET to MT torque
                Qratio=DqBET/DqMT;
                % Comparison of axial induction factors
                anew=0.5*(a+Tratio);
                % Comparison of tangential induction factors
                bnew=0.5*(b+Qratio);
                % Check for satisfaction of induction factor errors
                if (abs(anew-a)<1.0e-5),
                    if (abs(bnew-b)<1.0e-5),
                        TerminateSEQ=1;
                    end;
                end;
                % Reset induction factors
                a=anew;
                b=bnew;
                count=count+1;
                if (count>1000),
                    TerminateSEQ=1;
                end;
            end;
            % Total thrust and torque
            T=T+DtBET;
            Q=Q+DqBET;
        end;
        Tis(i,k)=T;
    end
end