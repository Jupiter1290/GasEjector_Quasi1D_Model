% This code designs the gas-gas ejector
%
%                     -  |
%      |-           -    |
%      |  -        -     |
%      |     - - -       |
%    1 |  2  |  5|       |6
%      |     - - -       |
%      |  -  | 4 |  -    |
%      | -   |   |    -  |
%            |___|      -|
%             3
%
close all;
clear;
clc;
P_1 = 12.5;                   % bar
T_1 = 350.0;                  % K
%V_1=100;                       %m/s
R_univ = 8.31434;         % m^3*Pa/(mol*K)
MW_nh3 = 17;                 % kg/kmol
MW_flue = 29;               % kg/kmol
mu_nh3 = 12e-6;             % Pa-s
mu_flue = 1e-5;             % Pa-s
gamma_nh3 = 1.33;
gamma_flue = 1.4;
D_1 = 4*0.0254;             % m
D_2=0.065;                 %m
D_3 = 1.5*0.0254;             % m
m_1 = 4.75;                    % kg/s
P_3 = 12.0;                     % bar
T_3 = 2800;                  % K
f = 1*0.015;                  %Darcy Friction factor
f_fanning=f/4;                 %Fanning friction factor
%R=8.4242;    %in (m^3 bar)/(Kg-K)
Cp_nh3=(R_univ/(MW_nh3*10^(-3)))*(gamma_nh3)/(gamma_nh3-1);        %in J/K-kg
Cp_flue=(R_univ/(MW_flue*10^(-3)))*(gamma_flue)/(gamma_flue-1);      %in J/K-kg
%Cv_nh3=Cp_nh3/gamma_nh3;      %in J/Kg-K 
%Cv_flue=Cp_flue/gamma_flue;   %in J/Kg-K
%Inlet conditions:
R_nh3=R_univ/0.017;
rho_1=P_1/(R_nh3*T_1)*10^5;
V_1=m_1/(rho_1*(pi/4)*D_1^2);
a_star_1=(gamma_nh3*(R_univ/(MW_nh3*0.001))*T_1)^(1/2);
machNo_1=V_1/a_star_1;

%For the converging section of the motive gas:
l_con=0.6988;  %in m

%for computing dA_dx;
syms 'x';
D=@(x) ((D_2-D_1)/l_con)*x+D_1;
A=@(x) (pi/4)*D(x)^2;
dA_dx=diff(A(x),x);
dA_dx_result=matlabFunction(dA_dx);

%computing the MachNo.
dM_dx=@(x,M,gamma)((M*(2+(gamma-1)*M^2))/2)*(((-1/A(x))*(dA_dx_result(x))+0.5*gamma*M^2*f/D(x))/(1-M^2));
odeFun=@(x,M)dM_dx(x,M,gamma_nh3);
sol=ode45(odeFun,[0,l_con],machNo_1);
[x,M]=ode45(odeFun,[0,l_con],machNo_1);

%Function of pressre along the converging nozzle:
P_con=@(gamma,x)P_1*((machNo_1/deval(sol,x))*((((2+(gamma-1)*machNo_1^2))/(2+(gamma-1)*(deval(sol,x))^2))^(1/2))*(A(0)/A(x)));
%isentropic realtion between temperatures (general relationship):
T=@(gamma,T_k,M,M_k)T_k*(2+(gamma-1)*M_k^2)/(2+(gamma-1)*M^2); %M->M^2
%Properties at point 2:
machNo_2=deval(sol,l_con);
P_2=P_con(gamma_nh3,l_con);
%Adiabatic relations:
T_2=T(gamma_nh3,T_1,deval(sol,l_con),machNo_1);
P_adiabatic=@(gamma,P_k,M,M_k)P_k*(M_k/M)*((2+(gamma-1)*(M_k)^2)/(2+(gamma-1)*(M)^2))^(1/2);
%Plotting (test)
figure()
%To show M=1 line;
xline=linspace(0,l_con,100);
plot(x,M);
title("Mach number along X in the Converging Section")
hold on;
plot(xline,(ones(size(xline))));
xlabel('Distance along the duct (in m)');
ylabel('MachNo');
ylim([0,1.5]);

%Duct parameters:
l_duct=0.7;
D_duct=0.0254*1.5;

%Solving the machNos for 3 and 4;
eqn = @(M) [
    (D_duct/(4*f_fanning))*((1 / gamma_flue) * ((1 / M(1)^2) - (1 / M(2)^2)) + ((gamma_flue + 1) / (gamma_flue * 2)) * log((M(1)^2 * (1 + ((gamma_flue - 1) / 2) * M(2)^2)) / (M(2)^2 * (1 + ((gamma_flue - 1) / 2) * M(1)^2))))-l_duct;
    (M(1) / M(2)) * (((2 + (gamma_flue - 1) * M(1)^2) / (2 + (gamma_flue - 1) * M(2)^2))^(1/2))-P_2/P_3;
];
guess=[1 1];
solutionMach=fsolve(eqn,guess);
machNo_3=solutionMach(1);
machNo_4=solutionMach(2);

%Mixing:
%Conservation of mass in the mixing chamber:
V_2=machNo_2*(gamma_nh3*(R_univ/(MW_nh3*10^(-3)))*T_2)^(1/2);
T_4=T(gamma_flue,T_3,machNo_4,machNo_3); %flag...
V_4=machNo_4*(gamma_flue*(R_univ/(MW_flue*10^(-3)))*T_4)^(1/2);
P_4_exp=P_adiabatic(gamma_flue,P_3,machNo_4,machNo_3); 
P_4=P_2;
rho_4=P_4/((R_univ/(MW_flue*10^(-3)))*T_4)*10^5;
m_4=rho_4*(pi/4)*D_duct^(2)*V_4;
n_flue=m_4/(MW_flue*10^(-3));   %In moles/sec;
n_nh3=m_1/(MW_nh3*10^(-3));     %In moles/sec;
chi=n_nh3/(n_flue+n_nh3);
m_5=m_1+m_4;
mass_frac_nh3=m_1/m_5; %versionCopy did not use massFraction
D_5=D_2; %Dia of conv exit=Dia of div entry
%V_5=((m_1*V_2)+(m_4*V_4))/m_5; %momentum conservation on the mixing chamber CV;
Cp_mix=((mass_frac_nh3*Cp_nh3)+((1-mass_frac_nh3)*Cp_flue));
totEnergy_2=m_1*(Cp_nh3*T_2+(V_2^(2))/2);
totEnergy_4=m_4*(Cp_flue*T_4+(V_4^(2))/2);
totEnergy_5=totEnergy_4+totEnergy_2;
%T_mix=(totEnergy_5-(m_5*V_5^(2))/2)/(m_5*Cp_mix); %correction1
MW_mix=(mass_frac_nh3*MW_nh3)+((1-mass_frac_nh3)*MW_flue);   %Careful! in gram/mol!! versionCopy used chi
gamma_mix=(mass_frac_nh3*gamma_nh3)+((1-mass_frac_nh3)*gamma_flue); %versionCopy used chi
eqn_2=@(t)[totEnergy_5-m_5*(Cp_mix*t(1)+(t(2)^2)/2);
    m_5-((P_2*10^5)/((R_univ/(0.001*MW_mix))*t(1)))*A(l_con)*t(2);
    ];
guess_2=[1 1];
options = optimoptions('fsolve', 'MaxFunctionEvaluations',200);
solution_5=fsolve(eqn_2,guess_2,options);
T_mix=solution_5(1);
V_5=solution_5(2);
sonic_5=(gamma_mix*(R_univ/(MW_mix*10^(-3)))*T_mix)^(1/2);
machNo_5=V_5/sonic_5;
%The diverging nozzle:
P_5=P_2;
D_6=0.1174;
l_div=1;   %in metres (sirs value, mine is 0.5)
T_5=T_mix;
T_5_total=T_5+V_5^2/(2*Cp_mix);
%for computing dA_dx_div;
syms 's';
D_div=@(s) ((D_6-D_5)/l_div)*s+D_5;
A_div=@(s) (pi/4)*(D_div(s))^2;
dA_ds_div=diff(A_div(s),s);
dA_ds_result_div=matlabFunction(dA_ds_div);

%computing the MachNo.
dM_ds_div=@(s,M,gamma)((M*(2+(gamma-1)*M^2))/2)*(((-1/A_div(s))*(dA_ds_result_div(s))+0.5*gamma*M^2*f/D_div(s))/(1-M^2));
odeFun_div=@(s,M)dM_ds_div(s,M,gamma_mix);
sol_div=ode45(odeFun_div,[0,l_div],machNo_5);
[s,M_div]=ode45(odeFun_div,[0,l_div],machNo_5);

%Function of pressre along the diverging nozzle:
P_div=@(gamma,s)P_5*((machNo_5/deval(sol_div,s))*((((2+(gamma-1)*machNo_5^2))/(2+(gamma-1)*(deval(sol_div,s))^2))^(1/2))*(A_div(0)/A_div(s)));

%At point 6:
machNo_6=deval(sol_div,l_div);
P_6=P_div(gamma_mix,l_div);
T_6=T(gamma_mix,T_5,machNo_6,machNo_5);

figure()
plot(s,M_div);
hold on;
title("Mach number along X in the Diverging Section")
xline_div=linspace(0,l_div,150);
plot(xline_div,ones(size(xline_div)));
xlabel('Distance along the duct (in m)');
ylabel('MachNo');
ylim([0,1.5]);
%Total_temps:
V_3=machNo_3*(gamma_flue*(R_univ/(MW_flue*0.001))*T_3)^(1/2);
T_3_total=T_3+V_3^2/(2*Cp_flue);
T_1_total=T_1+V_1^2/(2*Cp_nh3);
V_6=machNo_6*(gamma_mix*(R_univ/(MW_mix*0.001))*T_6)^(1/2);
T_6_total=T_6+V_6^2/(2*Cp_mix);

%Total pressures:
P_3_total=((T_3_total/T_3)^(gamma_flue/(gamma_flue-1)))*P_3;
P_1_total=((T_1_total/T_1)^(gamma_nh3/(gamma_nh3-1)))*P_1;
P_6_total=((T_6_total/T_6)^(gamma_mix/(gamma_mix-1)))*P_6;

%Outlet Pipe
%syms machNo_7
l_out=D_6*7;
out_function=@(m) (1/gamma_mix)*(1/machNo_6^2-1/m^2)+((gamma_mix+1)/(2*gamma_mix))*log((machNo_6^2*(1+((gamma_mix-1)/2)*m^2))/(m^2*(1+((gamma_mix-1)/2)*machNo_6^2)))-(4*f_fanning*l_out)/D_6;
machNo_7=fzero(out_function,0.05);
T_7=T_6*((2+(gamma_mix-1)*machNo_6^2)/(2+(gamma_mix-1)*machNo_6^2));
V_7=machNo_7*(gamma_mix*(R_univ/(MW_mix*0.001))*T_7)^(1/2);
T_7_total=T_7+V_7^2/(2*Cp_mix);
P_7=P_6*machNo_6/machNo_7*(T_7/T_6)^(1/2);
P_7_total=((T_7_total/T_7)^(gamma_mix/(gamma_mix-1)))*P_7;
%Plots for Analysis:
x_sampling_div=linspace(0,(l_div),150);
x_sampling_con=linspace(0,(l_con),150);
machNo_div=deval(sol_div,x_sampling_div);
machNo_con=deval(sol,x_sampling_con);
Pressure_div = arrayfun(@(s) P_div(gamma_mix, s), x_sampling_div);
Pressure_con=arrayfun(@(x) P_con(gamma_nh3,x),x_sampling_con);
T_con=@(gamma,x) T(gamma,T_1,deval(sol,x),machNo_1);
T_div=@(gamma,x) T(gamma,T_5,deval(sol_div,x),machNo_5);
Temperature_con=arrayfun(@(x) T_con(gamma_nh3,x),x_sampling_con);
Temperature_div=arrayfun(@(x) T_div(gamma_mix,x),x_sampling_div);
%Temperature_con=T_con(gamma_nh3,x_sampling_con);
%Temperature_div=T_div(gamma_mix,x_sampling_div);

figure()
plot(x_sampling_con,Temperature_con);
title("Temperature vs X in Converging Section")
xlabel('Along X (in metres)');
ylabel("Temperature (in Kelvin)");
figure()
plot(x_sampling_div,Temperature_div);
title("Temperature vs X in Diverging Section")
xlabel('Along X (in metres)');
ylabel("Temperature (in Kelvin)");
figure()
plot(x_sampling_con,Pressure_con);
title("Pressure vs X in Converging Section")
xlabel('Along X (in metres)');
ylabel("Pressure (in bar)")
figure()
plot(x_sampling_div,Pressure_div);
title("Pressure vs X in Diverging Section") 
xlabel('Along X (in metres)');
ylabel("Pressure (in bar)")