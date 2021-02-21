% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%Short program to evaluate the accuracy of Stirling Approximation. The    %
%Stirling Approximation is used to estimate the factorial of very large   %
%numbers, and is given by ln(n!)=n*ln(n)-n+(1/2)*ln(2*pi*n). Then, the    %
%factorial is estimated by n!=(2*pi*n)^(1/2)*(n/exp)^n. The program       %
%computes the ratio R=n!/((2*pi*n)^(1/2)*(n/exp)^n) for n=10, n=100,      %
%n=1000, n=10000 to determine the accuracy of the approximation. As       %                                                                        
%expected, the formula approaches exact value of R as n increases.        %                                                                 
%Programmed by Anthony Knighton on 2/07/2021.                             %                                            
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
n=10;
nmax=10000;
while (n<nmax+1)
   stirling = n*log(n) - n + 0.5*log(2*pi*n);  
   exact = sum(log(1:n));  
   logR = stirling - exact
   R=exp(logR)
   fprintf('n= %8d, R= %20.15e\n',n,R)
   n=n*10
end