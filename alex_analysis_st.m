function [time,bprofile,dFprofile,avg,filtprofile1]=alex_analysis_st(profile)
%set number of frames;length of light in ms, time of event in ms, 
%length of stimulation in ms, number of points to skip for fit 
s0=400;
tlight=1000;
tevent=550;
pulse=150;
skip=16;
%separate profile into sweeps
nsweeps=length(profile)/s0;
sprofile0=zeros(s0,nsweeps);
for p=1:nsweeps
sprofile0 (:,p) = profile((1+(p-1)*s0):(p*s0));
end
%select time point of onset of light using brush tool
hfig=figure('Name','Select time point of onset of light');
p2=plot(sprofile0(:,1),'b');
xlim([0 100])
hbrush=brush;
set(hbrush,'Enable','on');
hstart = uicontrol('Position',[250 10 80 20],'String','Ok',...
              'Callback','uiresume(gcbf)');
uiwait(gcf); 
brushed_locs = get(p2, 'BrushData');
close(hfig);
onsetpos=find(brushed_locs==1);
tstart=onsetpos+skip;
%get end of light position
hfig=figure('Name','Select time point of end of light');
p2=plot(sprofile0(:,1),'b');
xlim([100 200])
hbrush=brush;
set(hbrush,'Enable','on');
hstart = uicontrol('Position',[250 10 80 20],'String','Ok',...
              'Callback','uiresume(gcbf)');
uiwait(gcf);
brushed_locs = get(p2, 'BrushData');
close(hfig);
endpos=find(brushed_locs==1);
tend0=endpos-2;
taq=tlight/(endpos-onsetpos);
eventpos=round(tevent/taq)-1-skip;
eventlength=round(pulse/taq)+5;

%correct for baseline
sprofile=zeros(s0,nsweeps);
for p=1:nsweeps
bckg = mean(sprofile0(1:10,p));
sprofile(:,p)=sprofile0(:,p)-bckg;
end
%get rid of time prev. to start and do average
sprofile(1:tstart,:)=[];
s=size(sprofile,1);
nsweeps=size(sprofile,2);
bprofile=zeros(size(sprofile));
dFprofile=zeros(size(sprofile));
%fit and bleach correct
time=[0*taq:taq:(taq*(s-1))];
tend=tend0-tstart;
for p=1:nsweeps
tempprofile=transpose(sprofile(:,p));
fitprofile1=[tempprofile(1:(eventpos-1)),tempprofile((eventpos+eventlength):tend)];
fitprofile1=transpose(fitprofile1);
fittime1=[time(1:(eventpos-1)),time((eventpos+eventlength):tend)];
fittime1=transpose(fittime1);
bleachfit1= fit(fittime1(:,1), fitprofile1,'exp2');
bleachfactor1=bleachfit1(time);
limits=[(sprofile(tend,p)-10) (sprofile(1,p)+5)];
fh=figure (p);
plot(time,sprofile(:,p),'b')
hold on
plot(time,bleachfactor1,'k')
xlim([0 tlight])
ylim([min(limits) max(limits)])
waitfor(fh)
correct=transpose(bleachfactor1);
bprofile(:,p)=sprofile(:,p)./bleachfactor1;
end

%%run until here if you need to delete sweeps from bprofile, then redefine nsweeps
%example:
%bprofile(:,21)=[];
%nsweeps=40;

%transform to dF/F
dFprofile=zeros(s,nsweeps);
for p=1:nsweeps
bline = bprofile((eventpos-10):eventpos,p);
F = mean (bline);
dFprofile(:,p) = (bprofile(:,p)-F)./F;
end

%plot all signel trials
figure
plot(time,dFprofile,'b')
%do average
avg=zeros(1,s);
for j=1:s
avg (j) = mean(dFprofile(j,:));
end
%averaging filter
coeff3 = ones(1, 3)/3;
delay = mean(grpdelay(coeff3,1)); 
filtprofile1 = filter(coeff3, 1, avg); 
filtprofile1(1:delay)=[];
filtprofile1(1)=filtprofile1(2);

%plot it all
realtimeofevent=eventpos*taq;
figure
hax=axes;
plot(time,avg,'b')
hold on
plot(time(1:(s-1)),filtprofile1,'r') 
SP=realtimeofevent;
xlim([0 tend*taq])
ylim([-0.002 0.002])
line([SP SP],get(hax,'YLim'),'Color',[0 0 0])

%save stuff
folder_name = uigetdir;
oldFolder = cd(folder_name);
csvwrite('time.csv',time);
csvwrite('realtimeofevent.csv',realtimeofevent);
csvwrite('avg.csv',avg);
csvwrite('bprofile.csv',bprofile);
csvwrite('dFprofile.csv',dFprofile);
csvwrite('filtbavg1.csv',filtprofile1);
%save trasposed version for IGOR
csvwrite('dF_avg_igor.csv',transpose(avg));
csvwrite('filtbavg1_igor.csv',transpose(filtprofile1));
csvwrite('filtdF_avg_igor.csv',transpose(avg));
cd(oldFolder); 

