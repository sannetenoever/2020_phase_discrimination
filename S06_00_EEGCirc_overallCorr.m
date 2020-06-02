clear
addpath('<fieldtriplocation>');
addpath(genpath('<circstatloc>'));
locB = '<preprocDataLoc>';
ppnames = {'pp1';'pp2';'pp3';'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13';'pp14';'pp15';'pp16';'pp17';'pp18';'pp19';'pp20';'pp21'};
ID = 'R07_00'; % ID of original script

load([locB '/R01_02_stat'], 'stSI','SI');

%% load behavioral data and EEGcircCor stat
foi = 2:0.1:15;
load([locB 'S01_00_Behdata'], 'm');
load([locB 'R01_02_stat'], 'stSI','SI');

%% extract phase modulation index (average of all datapoints within the cluster):
clear v v2
for p = 1:21
    for fr = 1:2      
        d = squeeze((SI.CC(p,:,:)-median(SI.CCperm(p,:,:,:),4))./median(SI.CCperm(p,:,:,:),4));
        v(p,fr) = mean(d(stSI.posclusterslabelmat==fr));
    end
end

%% correlation
clear cors pv
close all
ppto = [1:21];
tnames = {'Categorizer theta';'Categorizer alpha'};
for fr = 1:2
    if fr == 1
        corv = v(ppto,2); % alpha
    elseif fr == 2
        corv = v(ppto,1); % theta
    end
    [cors(fr) pv(fr)] = corr(corv, m.slope(ppto)', 'tail', 'right');
    subplot(1,2,fr);
    lim = [min(corv) max(corv)];
    lim = [lim(1)-0.1 lim(2)+0.1];
    
    scatter(corv, m.slope(ppto)', 200, '.');
    hold on
    Fit = polyfit(corv,m.slope(ppto)',1);
    y = polyval(Fit,lim);
    plot([lim(1)-0.1 lim(2)+0.1],y);
    set(gca, 'xlim', [lim(1)-0.2 lim(2)+0.2]);
    xlabel('phase modulation');ylabel('slope');
    title(tnames{fr});
    text(lim(1)-0.2+0.08*[lim(2)+0.2-lim(1)-0.2], 0.039, ['r = ' num2str(round(cors(fr)*100)./100)]);
    text(lim(1)-0.2+0.08*[lim(2)+0.2-lim(1)-0.2], 0.0355, ['p = ' num2str(round(pv(fr)*100)./100)]);
end
set(gcf, 'position', [120 650 460 300]);


