% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%Program to numerically evaluate the first order derivative of sin(x)     %
%using the mean finite difference method. The numerical result is compared%
%with the exact result.                                                   %
%                                                                         %
%Programmed by Anthony Knighton on 02/04/2021                             %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
f = @(x) sin(x)
tollerance = 0.00001
N = 101;
dx = 2*pi/(N-1);
x = (0:N-1)*dx;
h = 2;
error = tollerance + 1;
fprintf('%5s %20s\n','h','maximum error')
while error>tollerance
    h = h/2;
    numerical = (sin(x+h)-sin(x-h))/(2*h);
    exact = cos(x);
    error = max(abs(numerical - exact));
    fprintf('%8.6f  %15.6e\n',h,error)
end
fprintf('Tollerance is met when h=%10f\n',h)
subplot(1,2,1)
p=plot(x,numerical,'o',x,exact);
set(p(1),'Linewidth',2,'Color','blue')
set(p(2),'Linewidth',2,'Color','green')
xlabel('x','Fontsize',14)
ylabel('$\frac{d}{dx} \sin(x)$','Interpreter','LaTeX','Fontsize',14) 
legend('numerical','exact')
legend('Location','northeast')
subplot(1,2,2)
q=plot(x,abs(numerical-exact));
set(q(1),'Linewidth',2,'Color','black')
xlabel('x','Fontsize',14)
ylabel('Absolute error','Fontsize',14)