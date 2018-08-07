x0=pi/4;
tol=10^(-6);
fp=abs(fprime(x0));
xvec=linspace(0,1,100);
fprimevec=-2*sin(xvec).*(2*cos(xvec)-2)+cos(xvec).*(sin(xvec)-1);
while fp>tol
    plot(xvec,fprimevec,'LineWidth',3,'Color','k')
    hold on;
    plot(xvec*0+x0,fprimevec,'LineWidth',2,'LineStyle',':','Color','k')
    hold on;
    plot(xvec,f2prime(x0)*(xvec-x0)+fprime(x0),'LineWidth',1,'Color','r')
    hold on;
    x0=x0-fprime(x0)/f2prime(x0);
    plot(xvec*0+x0,fprimevec,'LineWidth',2,'LineStyle',':','Color','b')
    hold off;
    x0=x0-fprime(x0)/f2prime(x0);
    fp=abs(fprime(x0));
    ylim([min(fprimevec,[],2),max(fprimevec,[],2)])
    legend('Derivada','x0','2da derivada','x1');
    grid on;
    pause;
end
display(['La solucion es x*=',num2str(x0)])
display(['La derivada en x*=',num2str(fprime(x0))])
display(['La 2da derivada en x*=',num2str(f2prime(x0))])
display(['La funcion en x*=',num2str((2*cos(x0)-2)^2+(sin(x0)-1)^2)])
function f=fprime(x)
f=-2*sin(x)*(2*cos(x)-2)+cos(x)*(sin(x)-1);
end
function f=f2prime(x)
f=-2*cos(x)*(2*cos(x)-2)+4*cos(x)*sin(x)...
    -sin(x)*(sin(x)-1)+cos(x)*(cos(x));
end
