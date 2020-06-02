clear
addpath('<fieldtriplocation>');
addpath(genpath('<circstatloc>'));
locB = '<preprocDataLoc>';
ppnames = {'pp1';'pp2';'pp3';'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13';'pp14';'pp15';'pp16';'pp17';'pp18';'pp19';'pp20';'pp21'};
ID = 'R01_03'; % ID of original script

load([locB '/R01_02_stat'], 'stSI','SI');

%% follow up statistics for maximum frequency/channel in cluster
clus = 2; % clus 1 = alpha, clus 2 = theta
Mask = stSI.stat;
Mask(stSI.posclusterslabelmat~=clus) = nan;
[m Finx] = max(nanmean(Mask,1)); % get maximum frequency in cluster
[m Cinx] = max(nanmean(Mask,2)); % get maximum channel in cluster
stSI.freq(Finx)

for pp = 1:length(ppnames)
    load([locB '/S03_00_phaseI_' ppnames{pp}], 'allf');
    % repeat for each sound type separately
    for it = 1:12
        inxA = find(allf.trialinfo(:,3)==1 & allf.trialinfo(:,1) == it);
        inxB = find(allf.trialinfo(:,3)==2 & allf.trialinfo(:,1) == it);
        I.ps.avg(pp,it,1) = circ_mean(angle(allf.fourierspctrm(inxA,Cinx,Finx)));
        I.ps.avg(pp,it,2) = circ_mean(angle(allf.fourierspctrm(inxB,Cinx,Finx)));  
        I.ps.avg(pp,it,3) = circ_dist(I.ps.avg(pp,it,1), I.ps.avg(pp,it,2));
        I.ps.avgI(pp,it) = length(inxA)./(length(inxA) + length(inxB)) > 0.2 & length(inxA)./(length(inxA) + length(inxB)) < 0.8;
    end    
    
    inxA = find(allf.trialinfo(:,3)==1);
    inxB = find(allf.trialinfo(:,3)==2);
    I.ind{pp,1} = angle(allf.fourierspctrm(inxA,Cinx,Finx));
    I.ind{pp,2} = angle(allf.fourierspctrm(inxB,Cinx,Finx));
    I.indT(pp,1) = circ_rtest(I.ind{pp,1}); % consistency of sound A
    I.indT(pp,2) = circ_rtest(I.ind{pp,2}); % consistency of sound B

    % Phase opposition test
    for p = 1:1000
       AV = angle(allf.fourierspctrm([inxA; inxB], Cinx,Finx));
       x = randperm(length(AV));
       PO(p) = circ_r(AV(x(1:length(inxA))))+circ_r(AV(x(length(inxA)+1:end)))-2*circ_r(AV);
    end
    POall = circ_r(AV(1:length(inxA)))+circ_r(AV(length(inxA)+1:end))-2*circ_r(AV);
    I.indT(pp,4) = length(find(PO>POall))/p;
    if I.indT(pp,4) == 0
        I.indT(pp,4) = 1/p;
    elseif I.indT(pp,4) == 1
        I.indT(pp,4) = (p-1)/p;
    end
    
    I.avg(pp,1) = circ_mean(I.ind{pp,1});
    I.avg(pp,2) = circ_mean(I.ind{pp,2});    
    I.avgdiff(pp) = circ_dist(I.avg(pp,1), I.avg(pp,2));
end
[I.test(1) z] = circ_rtest(I.avg(:,1)); % consistency of sound A
I.test(2) = circ_rtest(I.avg(:,2)); % consistency of sound B
[I.test(3) z] = circ_rtest(I.avgdiff); % consistency of difference
zv = norminv(I.indT(:,4)); % extract for p-value of phase opposition
[h p ci st] = ztest(zv,0,1)
I.test(4) = p;

%% circular plots for an individual
close all
pp = 7;
figure
subplot(1,2,1)
circ_plot(I.ind{pp,1},'hist',[],15,true,true,'linewidth',2,'color','r');
subplot(1,2,2)
circ_plot(I.ind{pp,2},'hist',[],15,true,true,'linewidth',2,'color','r');

%% circular plot over participants
ab = 15;
figure
subplot(1,3,1)
circ_plot(I.avg(:,1),'hist',[],ab,true,true,'linewidth',2,'color','r');
title('sound A')
subplot(1,3,2)
circ_plot(I.avg(:,2),'hist',[],ab,true,true,'linewidth',2,'color','r');
title('sound B')
subplot(1,3,3)
circ_plot(I.avgdiff','hist',[],ab,true,true,'linewidth',2,'color','r');
title('sound B-A')

%% plot separately for the different sound types
ab = 10;
figure
for it = 1:12
    subplot(3,4,it)
    inx = I.ps.avgI(:,it);
    circ_plot(I.ps.avg(inx,it,3),'hist',[],ab,true,true,'linewidth',2,'color','r');
    title(['sound B-A ' num2str(sum(inx))])
end
