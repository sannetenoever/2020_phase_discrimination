clear
addpath('<fieldtriplocation>');
addpath(genpath('<circstatloc>'));
locB = '<preprocDataLoc>';
ppnames = {'pp1';'pp2';'pp3';'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13';'pp14';'pp15';'pp16';'pp17';'pp18';'pp19';'pp20';'pp21'};
ID = 'R01_02'; % ID of the original script

%%
% load behavioral data
load([locB 'S01_00_Behdata'], 'm');
JNDold = 21.76;
JND = 2/3*JNDold;
freqtrA = logspace(log10(200-1.8*JND), log10(200-0.5*JND),6);
freqtrB = logspace(log10(200+0.5*JND), log10(200+1.8*JND),6);
freqtrTotal = [freqtrA freqtrB];

w = warning('query','last');
warning('off',w.identifier);

nperm = 1000;
for pp = 1:length(ppnames)
    
    %% load data and preallocate
    load([locB '/S03_00_phaseI_' ppnames{pp}], 'allf');
    
    % some preallocation
    SI.CC = zeros(length(allf.freq), length(allf.label));
    SI.CCperm = zeros(length(allf.freq), length(allf.label), 1000);
    ACC.CC = SI.CC;
    ACC.CCperm = SI.CCperm;
    
    %% extract the accuracies
    allf.trialinfo(:,7) = 0;
    inx = (allf.trialinfo(:,3) == 1 & allf.trialinfo(:,1) <= 6) | (allf.trialinfo(:,3) == 2 & allf.trialinfo(:,1) >= 6); 
    allf.trialinfo(inx,7) = 1;
    
    clear X
    RP = bootstrp(1000,@(x) x, 1:size(allf.trialinfo,1));
    Y = allf.trialinfo(:,3)-1;
    Yacc = allf.trialinfo(:,7);
    
    %% GLM to extract residuals of stimulus influence
    stimin = unique(allf.trialinfo(:,1));
    for st = 1:length(stimin)-1 % amount of dummy variables -1 
        X(:,st) = double(allf.trialinfo(:,2)==stimin(st));
        X(X(:,st)==0,st) = -1;
    end
    if length(stimin) > 1
        mdl = fitglm(X,Y,'Distribution','binomial');
        res = mdl.Residuals.Raw;
        % for accuracy
        mdl = fitglm(X,Yacc,'Distribution','binomial');
        resA = mdl.Residuals.Raw;
    else
       res = Y; 
       resA = Yacc;
    end
    
    %% circular-linear correlations       
    % for sound identity:
    SI.CC = circ_corrcl_mat(angle(allf.fourierspctrm), res);
    T = arrayfun(@(x) circ_corrcl_mat(angle(allf.fourierspctrm), res(randperm(size(allf.trialinfo,1)))), 1:nperm, 'uniformoutput', 0);
    T = cell2mat(T);
    SI.CCperm = reshape(T, [28 length(allf.freq) nperm]);       
      
    % for accuracy
    ACC.CC = circ_corrcl_mat(angle(allf.fourierspctrm), resA);
    T = arrayfun(@(x) circ_corrcl_mat(angle(allf.fourierspctrm), resA(randperm(size(allf.trialinfo,1)))), 1:nperm, 'uniformoutput', 0);
    T = cell2mat(T);
    ACC.CCperm = reshape(T, [28 length(allf.freq) nperm]);

    SI.label = allf.label;
    SI.freq = allf.freq;
    ACC.label = allf.label;
    ACC.freq = allf.freq;
    %save([locB '/' ID ppnames{pp} '_CC'], 'SI', 'ACC');
end

%% load and put into a single file
clear SI ACC
for pp = 1:length(ppnames)   
    sp = load([locB '/' ID ppnames{pp} '_CC'], 'SI', 'ACC');
    if pp == 1
       SI.label = sp.SI.label;
       SI.freq = sp.SI.freq;
       ACC= SI;
    end
    SI.CC(pp,:,:) = sp.SI.CC;
    SI.CCperm(pp,:,:,:) = sp.SI.CCperm;
    ACC.CC(pp,:,:) = sp.ACC.CC;
    ACC.CCperm(pp,:,:,:) = sp.ACC.CCperm;
end

%% statistics
load([locB '/S03_00_phaseI_' ppnames{1}], 'allf');
SI.dimord = 'subj_chan_freq';
SI.label{16} = 'Fp1';
SI.label{17} = 'Fp2';

V1 = SI;
V1.powspctrm = SI.CC;
V2 = SI;
V2.powspctrm = squeeze(median(SI.CCperm,4));
V3 = SI;
V3.powspctrm = ACC.CC;
V4 = SI;
V4.powspctrm = squeeze(median(ACC.CCperm,4));

cfg = [];
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.method = 'montecarlo';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;

cfgn = [];
cfgn.method = 'distance';
cfgn.neighbourdist = 0.25;
cfg.neighbours = ft_prepare_neighbours(cfgn, allf);
cfg.tail = 1;
cfg.clustertail = 1;
cfg.alpha = 0.05;
cfg.numrandomization = 1000;
cfg.design(1,1:2*length(ppnames))  = [ones(1,length(ppnames)) 2*ones(1,length(ppnames))];
cfg.design(2,1:2*length(ppnames))  = [1:length(ppnames) 1:length(ppnames)];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stSI = ft_freqstatistics(cfg, V1,V2); % response choice
stACC = ft_freqstatistics(cfg, V3,V4); % accuracy
%save([locB '/' ID '_stat'], 'stSI','stACC', 'SI','ACC');

%% plot
close all
load([locB '/' ID '_stat'], 'stSI','stACC', 'SI','ACC');

pos = [500 315 260 150];
pos2 = [500 315 200 150];

for clus = 1:2    
    load([base '/Programs/Matlab/fieldtrip-20120904/fieldtrip-20161231/template/layout/easycapM1']);
    stSI.CorrDiff = permute(squeeze(mean(SI.CC-median(SI.CCperm,4),1)), [2 1])';
    stSI.map = double(stSI.posclusterslabelmat==clus);
    stACC.CorrDiff = permute(squeeze(mean(ACC.CC-median(ACC.CCperm,4),1)), [2 1])';
    stACC.map = double(stACC.posclusterslabelmat==clus);
   
    figure
    inC = find(sum(stSI.map,2) > 0);
    inF = find(sum(stSI.map,1) > 0);
    CorrDiff = squeeze(mean(SI.CC(:,inC,:)-median(SI.CCperm(:,inC,:,:),4),1));    
    shadedErrorBar(stSI.freq, mean(CorrDiff,1),std(CorrDiff,[],1)./sqrt(size(CorrDiff,1)));
    hold on
    set(gca, 'xlim', SI.freq([1 end]), 'ylim', [-0.01 0.03]);
    ylim = get(gca,'ylim');
    x = area(SI.freq(inF), ones(1,size(inF,2)).*ylim(2));
    x.FaceAlpha = 0.1; x.EdgeAlpha = 0; x.FaceColor = [1 0 0];
    x = area(SI.freq(inF), ones(1,size(inF,2)).*ylim(1));
    x.FaceAlpha = 0.1; x.EdgeAlpha = 0; x.FaceColor = [1 0 0];
    shadedErrorBar(stSI.freq, mean(CorrDiff,1),std(CorrDiff,[],1)./sqrt(size(CorrDiff,1)));
    set(gcf, 'position',pos2);
    
    % topo of relevant interval
    figure
    cfg = [];
    cfg.highlight ='yes';
    cfg.highlightchannel = inC;
    cfg.xlim = stSI.freq(inF([1 end]));
    cfg.parameter = 'CorrDiff';
    cfg.layout = lay;
    cfg.interactive = 'yes';
    ft_topoplotER(cfg,stSI);
    colorbar
    set(gca, 'clim', [0 0.03]);
    
    set(gcf, 'position',pos);    
end


