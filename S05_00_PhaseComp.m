clear
addpath('<fieldtriplocation>');
addpath(genpath('<circstatloc>'));
locB = '<preprocDataLoc>';
ppnames = {'pp1';'pp2';'pp3';'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13';'pp14';'pp15';'pp16';'pp17';'pp18';'pp19';'pp20';'pp21'};
ID = 'R06_00'; % ID of original script

load([locB '/R01_02_stat'], 'stSI','SI');

%% extract the phases
clus = 2;
Mask = stSI.stat;
Mask(stSI.posclusterslabelmat~=clus) = nan;
[m Finx] = max(nanmean(Mask,1)); % maximum frequency
[m Cinx] = max(nanmean(Mask,2)); % maximum channel
stSI.freq(Finx);

load([locB '/S04_01_stat_freqbased'], 'stat', 'V1', 'V2');
Mask = stat{clus}.stat;
Mask(stat{clus}.posclusterslabelmat > 0);
[tp inx] = max(nanmean(Mask(Cinx,:),1)); % maximum time point

for pp = 1:length(ppnames)
    load([locB '/S03_00_phaseI_' ppnames{pp}], 'allf');    
    % phase for behavior
    inxA = find(allf.trialinfo(:,3)==1); % sound A
    inxB = find(allf.trialinfo(:,3)==2); % sound B
    I.Beh(pp,1) = circ_mean(angle(allf.fourierspctrm(inxA,Cinx,Finx))); % mean phase A
    I.Beh(pp,2) = circ_mean(angle(allf.fourierspctrm(inxB,Cinx,Finx))); % mean phase B
    I.Behdiff(pp) = circ_dist(I.Beh(pp,1), I.Beh(pp,2)); % phase difference A + B
    
    % phase for TRF    
    load([locB 'S03_00_phaseI_' ppnames{pp}], 'allf');
    load([locB 'S04_01_mTRF_' ppnames{pp}], 'mTRF');
    tps = 1:5:size(mTRF,3);
    mTRF = mTRF(:,:,tps,:);    
    mTRdiff = squeeze(mTRF(:,2,inx,Cinx)-mTRF(:,1,inx,Cinx));
    
    % median split mTRdiff
    inxA = find(mTRdiff < median(mTRdiff));
    inxB = find(mTRdiff > median(mTRdiff));
    I.TRF(pp,1) = circ_mean(angle(allf.fourierspctrm(inxA,Cinx,Finx)));
    I.TRF(pp,2) = circ_mean(angle(allf.fourierspctrm(inxB,Cinx,Finx)));
    I.TRFdiff(pp) = circ_dist(I.TRF(pp,1), I.TRF(pp,2));
            
    I.com(pp,1) = circ_dist(I.Beh(pp,1), I.TRF(pp,1)); % distance sound A
    I.com(pp,2) = circ_dist(I.Beh(pp,2), I.TRF(pp,2)); % distance sound B
end

%%
close all
ab = 15;
figure

circ_plot(I.com(:),'hist',[],ab,true,true,'linewidth',2,'color','r');       
[p z] = circ_vtest(I.com(:), 0);

