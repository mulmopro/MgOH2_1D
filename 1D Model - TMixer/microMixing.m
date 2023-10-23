function [time,varianceDissipation] = microMixing(s)

% this refers to T-mixer length of 0.04 m
x_=[0.00125	0.0025	0.00375	0.005	0.00625	0.0075	0.00875	0.01	0.01125	0.0125	0.01375	0.015	0.01625	0.0175	0.01875	0.02	0.02125	0.0225	0.02375	0.025	0.02625	0.0275	0.02875	0.03	0.03125	0.0325	0.03375	0.035	0.03625	0.0375	0.03875];
eps_=zeros(1,length(x_));
k_=[20.8744	17.7926	13.9368	10.8159	8.64852	7.04866	5.92258	5.08647	4.49663	3.97668	3.55206	3.21123	2.94116	2.72022	2.55403	2.39862	2.26849	2.15877	2.06601	1.98733	1.92699	1.86965	1.82107	1.77992	1.7451	1.71567	1.69309	1.67182	1.65394	1.63879	1.6262];

time = x_/s.velocity;
scaleTime=time(1);
time=time-scaleTime;

for i=1:length(time)
    eps_(i)=epsilon(time(i),s);
end

varianceDissipation=-s.Cphi./2.*eps_./k_;

% %   Plotting k and epsilon for checking the behaviour
% 
%     figure(1)
%     yyaxis left
%     plot(time,eps_,'LineWidth',2)
%     xlabel('time, s',"FontSize",25)
%     ylabel('$epsilon,$\ $\frac{m^2}{s^3}$',"FontSize",25,'Interpreter',"latex")
%     yyaxis right
%     plot(time,k_,'LineWidth',2)
%     ylabel('$k,$\ $\frac{m^2}{s^2}$',"FontSize",25,'Interpreter',"latex")
%     title('Micromixing',"FontSize",25)
% 
%     figure(2)
%     plot(time,-varianceDissipation)
%     xlabel('time, s')
%     ylabel('variance dissipation')
end