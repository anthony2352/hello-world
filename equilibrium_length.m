% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Problem 8.2: Program to find length of springs at mechanical            %         
% equlilbrium by solving system of linear equations using linsolve        %                                                                                
%                                                                         %
% Programmed by Anthony Knighton on 2/24/2021                             %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Define parameters
d=8;
l1=1;
l2=2;
l3=1;
l4=2;
k1=2;
k2=4;
k3=4;
k4=2;
% System of Equations
% dU/dx1=k1*(x1-l1)-k2*(x2-x1-l2)=0
% dU/dx2=k2*(x2-x1-l2)-k3*(x3-x2-l3)=0
% dU/dx3=k3*(x3-x2-l3)-k4*(d-x3-l4)=0
% Define coefficient array
A = [ k1+k2 -k2 0; -k2 k2+k3 -k3; 0 -k3 k3+k4];
% Define constant array
B = [ k1*l1-k2*l2; k2*l2-k3*l3; k3*l3+k4*d-k4*l4];
% Calculate Solution array
X = linsolve(A,B);
% Show the results
fprintf('x1 =%.6f\n',X(1))
fprintf('x2 =%.6f\n',X(2))
fprintf('x3 =%.6f\n',X(3))