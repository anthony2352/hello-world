% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%Program to evaluate the roots of a*x^2+x+1/4=0, where a(n+1)=a(n)/10.    %
%After a certain number of iterations, "a" hits machine epsilon. The exact%
%solution for a=0 is x=-.25. As "a" converges to machine epsilon, equation%
%x1 fails because eventually, b^2>>a*c, which leads to catastrophic       %
%cancellation. The error becomes particularly severe when a<<b because the%
%denominator becomes very small and the calculation approaches 0/0. This  %
%situation is mitigated by x2=-2*c/(b+(b^2-4*a*c)^(1/2)). Observe that x1 %
%eventually fails past a certain threshold.                               %
%                                                                         %
%Programmed by Anthony Knighton on 02/04/2021                             %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear 
b=1; 
c=0.25;
n=1;
a(1)= 1
while (a>eps())
  d = b^2-4*a(n)*c;
  x1(n) = (-b+sqrt(d))/(2*a(n));
  x2(n) = -2*c/(b+sqrt(d));
  fprintf('a= %22.16e, lhs=%22.16e, rhs=%22.16e \n',a(n),x1(n),x2(n));
  n=n+1;
  a(n)=a(n-1)/10;
end
%plot of root evaluated by original equation and plot of root evaluated by
%equation suitable for numerical calculation.
p=loglog(a(1:n-1),abs(x2+0.25),'o',a(1:n-1),abs(x1+0.25),'s');
set(p(1),'Linewidth',2,'Color','red')
set(p(2),'Linewidth',2,'Color','blue')
xlabel('a','Fontsize',14)
ylabel(texlabel('|x - x_0|'),'Fontsize',14)  % x_0 = root at a=0
legend('lhs','rhs')
legend('Location','Southeast')