% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%Program to determine the trajectory of a spring pendulum system.         %
%                                                                         %
%Programmed by Anthony Knighton on 02/21/2021                             %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear all;
% input parameter
a = input('intput value of a =>');
% independent variables
% q1 = L
% q2 = dL/dt
% q3 = theta
% q4 = dtheta/dt = omega
% All ODEs.
Q1 = @(Q1,Q2,Q3,Q4) Q2;
Q2 = @(Q1,Q2,Q3,Q4) Q1*Q4^2 - (Q1-1) + a*cos(Q3);
Q3 = @(Q1,Q2,Q3,Q4) Q4;
Q4 = @(Q1,Q2,Q3,Q4) -Q2*Q4/Q1 - a/Q1*sin(Q3);
% initial conditions
q1(1)=1.0;
q2(1)=0.0;
q3(1)=0.1;
q4(1)=0.0;
t(1)=0;
% control parameter for Runge-Kutta Integration
tmax=200;
N=10000;
h=tmax/N;
% 4th order Runge-Kutta method to solve the ODEs
for n=1:N-1
    % Euler step
    k11=Q1(q1(n),q2(n),q3(n),q4(n));
    k21=Q2(q1(n),q2(n),q3(n),q4(n));
    k31=Q3(q1(n),q2(n),q3(n),q4(n));
    k41=Q4(q1(n),q2(n),q3(n),q4(n));
    % 2nd order Runge-Kutta method
    q1_mid = q1(n)+k11*h/2;
    q2_mid = q2(n)+k21*h/2;
    q3_mid = q3(n)+k31*h/2;
    q4_mid = q4(n)+k41*h/2;
    k12=Q1(q1_mid,q2_mid,q3_mid,q4_mid);
    k22=Q2(q1_mid,q2_mid,q3_mid,q4_mid);
    k32=Q3(q1_mid,q2_mid,q3_mid,q4_mid);
    k42=Q4(q1_mid,q2_mid,q3_mid,q4_mid);
    % Predictor-corrector step
    q1_mid = q1(n)+k12*h/2;
    q2_mid = q2(n)+k22*h/2;
    q3_mid = q3(n)+k32*h/2;
    q4_mid = q4(n)+k42*h/2;
    k13=Q1(q1_mid,q2_mid,q3_mid,q4_mid);
    k23=Q2(q1_mid,q2_mid,q3_mid,q4_mid);
    k33=Q3(q1_mid,q2_mid,q3_mid,q4_mid);
    k43=Q4(q1_mid,q2_mid,q3_mid,q4_mid);
    % 4th order Runge-Kutta method
    q1_end = q1(n)+k13*h;
    q2_end = q2(n)+k23*h;
    q3_end = q3(n)+k33*h;
    q4_end = q4(n)+k43*h;
    k14=Q1(q1_end,q2_end,q3_end,q4_end);
    k24=Q2(q1_end,q2_end,q3_end,q4_end);
    k34=Q3(q1_end,q2_end,q3_end,q4_end);
    k44=Q4(q1_end,q2_end,q3_end,q4_end);
    % update point for next iteration
    q1(n+1)=q1(n)+(k11+2*(k12+k13)+k14)*h/6;
    q2(n+1)=q2(n)+(k21+2*(k22+k23)+k24)*h/6;
    q3(n+1)=q3(n)+(k31+2*(k32+k33)+k34)*h/6;
    q4(n+1)=q4(n)+(k41+2*(k42+k43)+k44)*h/6;
    % update time parameter for next iteration
    t(n+1)=t(1)+n*h;
end
% plot angular coordinate (theta)
subplot(1,3,1);
q=plot(t,q3);
xlabel('t');
ylabel('theta');
set(q(1),'Color','red','Linewidth',2);
% plot angular coordinate (L)
subplot(1,3,2);
p=plot(t,q1);
xlabel('t');
ylabel(texlabel('length'));
set(p(1),'Color','red','Linewidth',2);
% trajectory of the bob in xy coordiates
x=q1.*sin(q3);
y=q1.*cos(q3)-1.0;
subplot(1,3,3);
plot(x,y);
xlabel('x');
ylabel('y');