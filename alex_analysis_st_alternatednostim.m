%% @Victoria Gonzalez-Sabater   -  2018

function [time,bprofile,dFprofile,avg,filtprofile1,factors]=alex_analysis_st_alternatednostim(profile)
%set number of frames;length of light in ms, time of event in ms, 
%length of stimulation in ms, number of points to skip for fit 
[s0 nsweeps]= size(profile);
tlight=3000;
tevent=550;
pulse=2000;
skip=40;
stimsweeps=[1:2:nsweeps];
nostimsweeps=[2:2:nsweeps];

parameters  = struct('tlight',tlight,'tevent',tevent,'pulse',pulse,'skip',skip)
%separate profile into sweeps

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
xlim([300 500])
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
bckg = mean(sprofile0(1:8,p));
sprofile(:,p)=sprofile0(:,p)-bckg;
end
%get rid of time prev. to start
sprofile(1:tstart,:)=[];
%separate stim and no stim
sprofile_nostim=sprofile(:,nostimsweeps);
sprofile_stim=sprofile(:,stimsweeps);
fh=figure;
subplot(1,2,1)
plot(sprofile_nostim)
subplot(1,2,2)
plot(sprofile_stim)
 promptMessage = sprintf('Do you want to Continue processing,\nor Cancel to abort processing?');
button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');
if strcmpi(button, 'Cancel')
	return
else
    close (fh)
end
tend=tend0-tstart;
%fit the no stim sweeps to obtain time constants
s=size(sprofile_nostim,1);
time=[0*taq:taq:(taq*(s-1))];
tempprofile=(mean(sprofile_nostim,2));
fitprofile1=tempprofile(1:tend);
fittime1=time(1:tend);
fittime1=transpose(fittime1);
bleachfit1= fit(fittime1(:,1), fitprofile1,'exp2');
bleachfactor1=bleachfit1(time);
coefficients=coeffvalues(bleachfit1);
checkprofile=tempprofile./bleachfactor1;
%plot to check
fh=figure;
subplot(1,2,1)
plot(time,tempprofile,'b')
hold on
plot(time,bleachfactor1,'k')
ylim([(tempprofile(tend)-30) (tempprofile(1)+30)])
xlim ([time(1) (time(tend)+50)])
subplot(1,2,2)
plot(time,checkprofile)
waitfor('fh')

%use time constants to fit each single run
[bprofile,factors]=bleachcorrect_alex(sprofile_stim,taq,tend,eventpos,eventlength,coefficients(2),coefficients(4));
promptMessage = sprintf('Do you want to Continue processing,\nor Cancel to abort processing?');
button = questdlg(promptMessage, 'Continue', 'Continue', 'Cancel', 'Continue');
if strcmpi(button, 'Cancel')
	return
end

s=size(sprofile_stim,1);
time=[0*taq:taq:(taq*(s-1))];
nsweeps=size(bprofile,2);
%transform to dF/F
dFprofile=zeros(s,nsweeps);
for p=1:nsweeps
bline = bprofile((eventpos-10):eventpos,p);
F = mean (bline);
dFprofile(:,p) = (bprofile(:,p)-F)./F;
end

%plot all single trials
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
xlim([0 3000])
ylim([-0.1 0.1])
line([SP SP],get(hax,'YLim'),'Color',[0 0 0])

%calculate raw df subtractions

dfsprofile = df_nonfit(profile) % if you want to see direct subtraction

%save stuff
folder_name = uigetdir;
oldFolder = cd(folder_name);
csvwrite('profile.csv',profile);
csvwrite('time.csv',time);
csvwrite('realtimeofevent.csv',realtimeofevent);
csvwrite('avg.csv',avg);
csvwrite('bprofile.csv',bprofile);
csvwrite('dFprofile.csv',dFprofile);
csvwrite('filtbavg1.csv',filtprofile1);
%save trasposed version for IGOR
csvwrite('dF_avg_igor.csv',transpose(avg));
csvwrite('filtbavg1_igor.csv',transpose(filtprofile1));
%save fitting properties
csvwrite('fitfactors.csv',factors);
csvwrite('nostimfit.csv',coefficients);
csvwrite('dfsprofile.csv',dfsprofile);
save('parameters','parameters');
cd(oldFolder); 

