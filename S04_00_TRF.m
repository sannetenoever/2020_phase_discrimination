clear
addpath('<fieldtriplocation>');
addpath(genpath('<circstatloc>'));
addpath('mTRFToolboxLocation>');
locB = '<preprocDataLoc>';
ppnames = {'pp1';'pp2';'pp3';'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13';'pp14';'pp15';'pp16';'pp17';'pp18';'pp19';'pp20';'pp21'};
cd([base '/Experiments/Collaborations/PhasePred/03_Scripts']);
ID = 'S04_01'; % ID of original script

%% now do the mTRF based on ambiguous stimuli 
load([locB 'S01_00_Behdata'], 'm');
load([locB 'S04_00_SndEnv'],'sndenvI');
fs=1000;
tp = [-0.1 0.5995];

JNDold = 21.76;
JND = 2/3*JNDold;
freqtrA = logspace(log10(200-1.8*JND), log10(200-0.5*JND),6);
freqtrB = logspace(log10(200+0.5*JND), log10(200+1.8*JND),6);
freqtrTotal = [freqtrA freqtrB];

for pp = 1:21
    clear mTRF
    load([locB 'S02_00_eyeblinkcorrectionRMT_' ppnames{pp}], 'datA_RM');
    stimtouse = find(freqtrTotal >= m.th20(pp) & freqtrTotal <= m.th80(pp));
    cfg = [];
    cfg.toilim = tp;
    cfg.trials = find(ismember(datA_RM.trialinfo(:,1), stimtouse));
    datA_RM = ft_redefinetrial(cfg, datA_RM);   
    % for each single trial do the mTRF
    for it = 1:length(datA_RM.trial)
       sttouse = squeeze(sndenvI(datA_RM.trialinfo(it,8), datA_RM.trialinfo(it,7), 1,:)); 
       mTRF(it,1,:,:) = mTRFtrain(sttouse, zscore(datA_RM.trial{it}',[],2),fs,1,-100, 600,1000);                            
       sttouse = squeeze(sndenvI(datA_RM.trialinfo(it,8), datA_RM.trialinfo(it,7), 12,:)); 
       mTRF(it,2,:,:) = mTRFtrain(sttouse, zscore(datA_RM.trial{it}',[],2),fs,1,-100, 600,1000);                            
    end
    %save([locB '\S04_01_mTRF_' ppnames{pp}], 'mTRF');
end

%% circular correlation TRF difference and phase
load([locB 'S04_01_mTRF_' ppnames{pp}], 'mTRF');
tps = 1:5:size(mTRF,3); % just for speeding up things

for pp = 1:21
    CC = zeros(28,131,length(tps));
    CCpermT = zeros(28,131,length(tps));
    CCpermT2 = zeros(28,131,length(tps));
    % now do the circular-linear correlation for all frequency and time
    % points
    load([locB 'S03_00_phaseI_' ppnames{pp}], 'allf');
    load([locB 'S04_01_mTRF_' ppnames{pp}], 'mTRF');
    mTRF = mTRF(:,:,tps,:);
    % make the perms
    % for the angle perms:
    RA = arrayfun(@(x) randperm(size(mTRF,1)), [1:500], 'uniformoutput', 0);
    % for the difference perms (random times -1):
    RD = arrayfun(@(x) ((random('Discrete Uniform',2,[size(mTRF,1),1])-1).*2)-1, [1:500], 'uniformoutput', 0);   
   
    for tp = 1:length(tps)
        for ch = 1:28
            mTRdiff = squeeze(mTRF(:,2,tp,ch)-mTRF(:,1,tp,ch));
            CC(ch,:,tp) = circ_corrcl_mat(squeeze(angle(allf.fourierspctrm(:,ch,:))), mTRdiff); 
            % permutations on the angle
            x = cellfun(@(x)circ_corrcl_mat(squeeze(angle(allf.fourierspctrm(:,ch,:))), mTRdiff(x))',RA, 'uniformoutput',0); 
            CCpermT(ch,:,tp) = median(cell2mat(x),2);
            % permutations on the difference
            mTRdiffP = cellfun(@(x) mTRdiff.*x, RD, 'uniformoutput', 0);
            x = cellfun(@(x) circ_corrcl_mat(squeeze(angle(allf.fourierspctrm(:,ch,:))), x)', mTRdiffP,'uniformoutput',0);
            CCpermT2(ch,:,tp) = median(cell2mat(x),2);            
        end
    end 
    %save([locB '\S04_01_CC_' ppnames{pp}], 'CC', 'CCpermT','CCpermT2', 'tps');
end

%% load relevant info for statistics
clear CC CCpermT CCpermT2 CC2
for pp = 1:21
    x = load([locB '/S04_01_CC_' ppnames{pp}], 'CC', 'CCpermT','CCpermT2');
    CC(pp,:,:,:) = x.CC;
    CCpermT(pp,:,:,:) = x.CCpermT;
    CCpermT2(pp,:,:,:) = x.CCpermT2;
end

load([base '/Programs/Matlab/fieldtrip-20120904/fieldtrip-20161231/template/layout/easycapM1']);
load([locB '/S03_00_phaseI_' ppnames{1}], 'allf');
allf.label{16} = 'Fp1';
allf.label{17} = 'Fp2';

% for statistics
fftr.label = allf.label;
fftr.time = linspace(-0.1, 0.6, size(CC,4));
fftr.freq = allf.freq;
fftr.powspctrm = squeeze(mean(CC,1))- squeeze(mean(CCpermT,1));
fftr.dimord = 'chan_freq_time';

% load the frequency bands of interest
load([locB 'S03_00_stat'], 'stSI');

% check the significant bands:
freq{1} = stSI.freq(find(sum(stSI.posclusterslabelmat == 1) > 1));
freq{2} = stSI.freq(find(sum(stSI.posclusterslabelmat == 2) > 1));

%% statistics
clear V1 V2 stat
for cltol = 1:2
    tp = 1:length(fftr.time);
    frinx = find(fftr.freq >= freq{cltol}(1) & fftr.freq <= freq{cltol}(end));
     
    V1{cltol} = fftr;
    V1{cltol}.dimord = 'subj_chan_time';
    V1{cltol}.individual = squeeze(mean(CC(:,:,frinx,tp),3));
    V1{cltol}.label{16} = 'Fp1';
    V1{cltol}.label{17} = 'Fp2';
    V1{cltol}.time = fftr.time(tp);
    V1{cltol} = rmfield(V1{cltol}, 'freq');
    V2{cltol} = V1{cltol};
    V2{cltol}.individual = squeeze(mean(CCpermT(:,:,frinx,tp),3));
    datA_RM.label = V1{cltol}.label;
    
    cfg = [];
    cfg.statistic = 'ft_statfun_depsamplesT';
    cfg.method = 'montecarlo';
    cfg.correctm = 'cluster';
    cfg.clusterthreshold = 'nonparametric_individual';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan = 1;
    
    cfgn = [];
    cfgn.method = 'distance';
    cfgn.neighbourdist = 0.25;
    cfg.neighbours = ft_prepare_neighbours(cfgn, datA_RM);
    cfg.tail = 1;
    cfg.clustertail = 1;
    cfg.alpha = 0.05;
    cfg.numrandomization = 1000;
    cfg.design(1,1:2*length(ppnames))  = [ones(1,length(ppnames)) 2*ones(1,length(ppnames))];
    cfg.design(2,1:2*length(ppnames))  = [1:length(ppnames) 1:length(ppnames)];
    cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
    cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
    
    stat{cltol} = ft_timelockstatistics(cfg, V1{cltol},V2{cltol});
end
%save([locB '\S04_01_stat_freqbased'], 'stat', 'V1', 'V2');

%% plot
load([locB 'S04_01_stat_freqbased'], 'stat', 'V1', 'V2');

for cltol = 1:2
    figure
    in = find([stat{cltol}.posclusters.prob] < 0.05);
    stat{cltol}.clus = ismember(stat{cltol}.posclusterslabelmat, in);
    stat{cltol}.dimord = 'chan_time';
    stat{cltol}.stat = stat{cltol}.stat;
    for inc = 1:length(in)
        sumst = sum(stat{cltol}.stat.*(stat{cltol}.posclusterslabelmat == inc),2);
        [t chinx] = sort(sumst,'descend');
        chinx(t==0) = [];
        tw = find(sum(stat{cltol}.posclusterslabelmat == inc,1)>0);
        dtp = smoothdata(squeeze(mean(V1{cltol}.individual(:,chinx(1:5),:),2)-mean(V2{cltol}.individual(:,chinx(1:5),:),2)),2,'movmean',7);
         
        subplot(4,2,inc+(inc-1))
        shadedErrorBar(stat{cltol}.time, squeeze(mean(dtp,1)),squeeze(std(dtp,1))./sqrt(size(dtp,1))); hold on
        set(gca, 'ylim', [-0.01 0.03], 'xlim', [-0.1 0.6]);
        if cltol == 1
            set(gca, 'ylim', [-0.01 0.06]);
        end
        ylim = get(gca, 'ylim');
        x = area(stat{cltol}.time(tw), ones(1,size(tw,2)).*ylim(2));
        x.FaceAlpha = 0.1; x.EdgeAlpha = 0; x.FaceColor = [1 0 0];
        
        % now the topo  
        subplot(4,2,inc+(inc-1)+1)
        stat{cltol}.dtp = squeeze(mean(V1{cltol}.individual,1)-mean(V2{cltol}.individual,1));
        cfg = [];
        cfg.parameter = 'dtp';
        cfg.xlim = stat{cltol}.time(tw([1 end]));
        cfg.highlight = 'on';
        cfg.highlightchannel = chinx(1:5);
        cfg.highlightsymbol = '*';
        cfg.layout = lay;
        cfg.zlim = [-0.01 0.025];
        if cltol == 1
           cfg.zlim = [-0.01 0.045];
        end
        cfg.colorbar =  'yes';
        ft_topoplotER(cfg, stat{cltol});
    end       
    set(gcf, 'position',[-1500 -270 1100 900]);    
end


