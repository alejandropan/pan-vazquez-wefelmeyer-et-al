function [bcprofile,factors]=bleachcorrect_alex(sprofile,taq,tend,eventpos,eventlength,tau1,tau2)
s=size(sprofile,1);
nsweeps=size(sprofile,2);
time=[0:taq:taq*(s-1)];
fittime1=[time(1:(eventpos-1)),time((eventpos+eventlength):tend)];
fittime1=transpose(fittime1);
fun = @(x,xdata)x(1)*exp(tau1*xdata)+x(2)*exp(tau2*xdata);
x0=[300,0];
factors=zeros(2,nsweeps);
bcprofile=zeros(s,nsweeps);
for p=1:nsweeps
tempprofile=transpose(sprofile(:,p));
fitprofile1=[tempprofile(1:(eventpos-1)),tempprofile((eventpos+eventlength):tend)];
ydata=transpose(fitprofile1);
x = lsqcurvefit(fun,x0,fittime1,ydata);
factors(:,p)=x;
fit1=@(xdata)x(1)*exp(tau1*xdata)+x(2)*exp(tau2*xdata);
bleachfactor1=fit1(time);
correct=transpose(bleachfactor1);
bcprofile(:,p)=sprofile(:,p)./correct;
limits=[(sprofile(tend,p)-10) (sprofile(1,p)+5)];
fh=figure(p);
subplot(1,2,1)
plot(time,sprofile(:,p),'b')
hold on
plot(time,bleachfactor1,'k')
ylim([min(limits) max(limits)])
xlim ([time(1) (time(tend)+50)])
subplot(1,2,2)
plot(time,bcprofile(:,p))
xlim ([time(1) (time(tend)+50)])
ylim ([(bcprofile(1,p)-3) (bcprofile(1,p)+3)])
waitfor(fh)
end