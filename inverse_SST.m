% Any use of this software must refer to the publication:
% Ryan Borowiecki, Vadim A. Kravchinsky, Mirko van der Baan, Roberto Henry Herrera, 2023. 
% The Synchrosqueezing Transform to evaluate paleoclimate cyclicity. Computers and Geosciences, in press.

clc; close all; clear all; clc
data=csvread('synth_data.csv');
st=data(:,1); t=data(:,2);
% extract header info
n = length(st);
dt = 0.05;
t=(1:n)*dt; % time vector

sm=smooth(st,0.6,'loess');

figure
subplot(2,1,1)
plot(t,st,t,sm)
st=st-sm; % Remove Trend
subplot(2,1,2)


st=smooth(st,0.035,'loess'); % High Frequency
plot(t,st)

%% Computation
variance = std(st)^2;
st = (st - mean(st))/sqrt(variance) ;

% Spline Interpolation
tu=linspace(0,max(t),1024);
fn_resamp=interp1(t,st,tu,'spline');


%% Windowed
dt = tu(2)-tu(1) ;

% Calculate coi:
coi=coi_calc(fn_resamp,dt);

% Caulculate SST:
[Tff,w]=wsst(fn_resamp, 1./dt);


% Inverse SST Individual Components:
freqrange1=[1.75,2.5];
xrec1=iwsst(Tff,w,freqrange1);

freqrange2=[0.85,1.2];
xrec2=iwsst(Tff,w,freqrange2);

freqrange3=[0.62,0.80];
xrec3=iwsst(Tff,w,freqrange3);

freqrange4=[0.35,0.55];
xrec4=iwsst(Tff,w,freqrange4);

xrec=xrec1+xrec2+xrec3+xrec4;

figure
set(gcf, 'Position',  [0, 0, 1850, 600])
subplot(922)
plot(tu,xrec1,'b'); xlim([0 12])
%title(['Mode 1 reproduces ' num2str(floor(mode_v(1))) '% of Amplitude'])
set(gca,'XTicklabels',[]); set(gca,'fontsize',12);

subplot(924)
plot(tu,xrec2,'m'); xlim([0 12])
mode2_v=max(xcorr(fn_resamp,xrec2))./1000;%title(['Mode 2 reproduces ' num2str(floor(mode_v(2))) '% of Amplitude'])
set(gca,'XTicklabels',[]); set(gca,'fontsize',12);

subplot(926)
plot(tu,xrec3,'g'); xlim([0 12])
%title(['Mode 3 reproduces ' num2str(floor(mode_v(3))) '% of Amplitude'])
set(gca,'XTicklabels',[]); set(gca,'fontsize',12);

subplot(928)
plot(tu,xrec4,'c'); xlim([0 12])
%title(['Mode 4 reproduces ' num2str(floor(mode_v(4))) '% of Amplitude'])
set(gca,'fontsize',12); xlabel('Time (ka)')

total_var=floor((max(xcorr(fn_resamp,xrec))./100).^2);

subplot(921)
plot(t,st,tu,xrec); title(['Reproduces ' num2str(total_var) '% of Amplitude'])
axis([0 12 -2.5 2.5]); set(gca,'XTicklabels',[]); set(gca,'fontsize',12);


% Noise Red Injection:
sigma = 0.25 ;%sqrt( var(clean)*10.^( -snrdb /10 ) );
realization=100;

for i=1:realization
    i;
    noise=rednoise(length(tu));
    stack=fn_resamp'+(sigma.*noise);
    [Tf,w]=wsst(stack, 1./dt);
    s_sgs(:,:,i)=Tf;
end
    
s_conf=sum(s_sgs,3)./realization;


% Plot the Results:
subplot(9,2,[3,5,7])
imagesc(tu, -log2(w), abs(Tff));
colormap (flipud(bone)); set(gca,'YDir','reverse');
ylabel('Period (kyr)'); set(gca,'XTicklabels',[]);
yt=get(gca,'YTick'); yticks=2.^yt'; yticklabels({num2str(yticks)});
axis([0 12 -2 2]); set(gca,'fontsize',12); caxis([0 .075])

hold on
area(tu,log2(coi),2); alpha(0.5);
hold on
contour(tu,-log2(w),abs(s_conf),[-99,0.01],'r');
hold on
% Inverse SST Search areas
plot(tu,-log2(freqrange1(1)).*ones(length(tu)),'b')
hold on

plot(tu,-log2(freqrange1(2)).*ones(length(tu)),'b')
hold on

plot(tu,-log2(freqrange2(1)).*ones(length(tu)),'m')
hold on
plot(tu,-log2(freqrange2(2)).*ones(length(tu)),'m')
hold on

plot(tu,-log2(freqrange3(1)).*ones(length(tu)),'g')
hold on
plot(tu,-log2(freqrange3(2)).*ones(length(tu)),'g')
hold on

plot(tu,-log2(freqrange4(1)).*ones(length(tu)),'c')
hold on
plot(tu,-log2(freqrange4(2)).*ones(length(tu)),'c')
hold off


%% Ridge Extraction

% Calculate coi
coi=coi_calc(fn_resamp,dt);

% Calculate SST:
[Tff,w] = wsst(fn_resamp, 1./dt);

% Ridges
[fridge,iridge] = wsstridge(Tff,12,w,'NumRidges',3);

% Inverse SST:
xrec = iwsst(Tff,iridge);



%% Plot the Results:
subplot(9,2,[13,15,17])
imagesc(tu, -log2(w), abs(Tff));
colormap (flipud(bone)); set(gca,'YDir','reverse');
xlabel('Time (ka)'); ylabel('Period (kyr)');
yt=get(gca,'YTick'); yticks=2.^yt'; yticklabels({num2str(yticks)});
axis([0 12 -2 2]); set(gca,'fontsize',12); caxis([0 .05])
hold on
area(tu,log2(coi),2); alpha(0.5);
hold on
plot(tu,-log2(fridge),'r--','linewidth',2);

subplot(9,2,12)
mode1=1;
plot(tu,xrec(:,mode1))
mode1_v=max(xcorr(fn_resamp,xrec(:,mode1)))./1000;
%title(['Mode 1 reproduces ' num2str(floor(mode1_v*100)) '% of Variability'])
set(gca,'fontsize',12); set(gca,'XTicklabels',[]);

subplot(9,2,14)
mode2=2;
plot(tu,xrec(:,mode2))
mode2_v=max(xcorr(fn_resamp,xrec(:,mode2)))./1000;
%title(['Mode 2 reproduces ' num2str(floor(mode2_v*100)) '% of Variability'])
set(gca,'fontsize',12); set(gca,'XTicklabels',[]); 

subplot(9,2,16)
mode3=3;
plot(tu,xrec(:,mode3))
mode3_v=max(xcorr(fn_resamp,xrec(:,mode3)))./1000;
%title(['Mode 3 reproduces ' num2str(floor(mode3_v*100)) '% of Variability'])
set(gca,'fontsize',12); xlabel('Time (ka)')


total_var=floor(100*(mode1_v.^2+mode2_v.^2+mode3_v.^2));


subplot(9,2,11)
plot(tu,fn_resamp,tu,sum(xrec,2)); axis([0 12 -2.5 2.5]); set(gca,'xticklabel',[])
title(['Reproduces ' num2str(total_var) '% of Amplitude'])
set(gca,'fontsize',12);
