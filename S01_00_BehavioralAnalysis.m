% This file performs the behavioral analysis based on the logfiles
loc =  '<locLogFiles>';
modelfreeloc = '<modelfreeloc>';

cd(loc);
ppnames = {'pp1';'pp2';'pp3'; 'pp4';'pp5';'pp6';'pp7';'pp8';'pp9';'pp10';'pp11';'pp12';'pp13';'pp14';'pp15';'pp16';'pp17';'pp18';'pp19';'pp20';'pp21'};

%%
close all
% estimae the frequencies that we used in the experiment:
JNDold = 21.76;
JND = 2/3*JNDold;
freqtrA = logspace(log10(200-1.8*JND), log10(200-0.5*JND),6);
freqtrB = logspace(log10(200+0.5*JND), log10(200+1.8*JND),6);
freqtrTotal = [freqtrA freqtrB]; 

for pp = 1:length(ppnames)
    %%
    clear newV vars
    figure;
    cntplot = 1;
    
    %% loop to read in the data from the logfile
    alld = dir([ppnames{pp} '_train*.txt']);
    for afps = 1:length(alld)
        csv_file = fopen(alld(afps).name,'rt');
        formatSpec = '%s';
        C = textscan(csv_file, formatSpec, 'Delimiter', ';'); % reads in the data as a single column
        st = 9; % line to start reading from the logfile
        e = 20; % line at to end reading from the lgofiles
        if pp >= 12 % at participant 12 we added some extra information at the top of the logifles
            st = 11;
            e = 22;
        end
        if pp >= 13 % at participant 12 we added a column to the logfiles
            st = 11;
            e = 23;
        end
        vd = [1 2 3 4 6 7 8 9 10 11]; % which variables we will need to read in
        v = 5; % these indexes are string
        ll = e-st; % this is the stepsize (dependent on the amount of columns)
        for i = 1:max(vd)
            vars(afps).(C{1}{st+i-1}(1:2)) = C{1}(e+i:ll+1:end);
            if ismember(i,vd)
                vars(afps).(C{1}{st+i-1}(1:2)) = cellfun(@(x) str2double(x), vars(afps).(C{1}{st+i-1}(1:2)));
            end
            if ismember(i,v)
                vars(afps).(C{1}{st+i-1}(1:2)) = cellfun(@(x) str2double(x(1:end-2)), vars(afps).(C{1}{st+i-1}(1:2)));
            end
        end
    end
    
    %% compare the separate files into one big one (only keep relevatn variable)
    newV = [];
    v = [3 6 8 1 2 4 7 9]; % which variables to keep
    for i = 1:length(v)
        newV.(C{1}{st+v(i)-1}(1:2)) = [];
    end
    for afps = 1:length(alld)
        l = length(vars(afps).(C{1}{st+i-1}(1:2)));
        for i = 1:length(v)
            newV.(C{1}{st+v(i)-1}(1:2))(end+1:end+l) = vars(afps).(C{1}{st+v(i)-1}(1:2));
        end
    end
    
    %% get behavioral data, overall slope
    for pit = 1:12
        inx = find(newV.Pi == pit);
        if mod(pp,2) == 0 | pp == 1
            pcor(pp,pit) = sum(newV.Re(inx) == 1)./length(inx);
            amtri(pp,pit) = sum(newV.Re(inx)== 1);
        else
            pcor(pp,pit) = sum(newV.Re(inx) == 2)./length(inx);
            amtri(pp,pit) = sum(newV.Re(inx)== 2);
        end
        amtr(pp,pit) = length(inx);
    end
    % fit psychometric curve:
    addpath(modelfreeloc);
    % cumulative gaussian distribution
    degpol = 1; % Degree of the polynomial
    guessing = 0; % guessing rate
    lapsing = 0; % lapsing rate
    b = binomfit_lims(amtri(pp,:)', amtr(pp,:)', freqtrTotal', degpol, 'probit', guessing, lapsing );
    
    % Plot the fitted curve and normal curve
    pfitnew = binomval_lims(b, freqtrTotal');
    m.or = pcor(pp,:);
    
    numxfit = 199; % Number of new points to be generated minus 1
    xfit = [min(freqtrTotal):(max(freqtrTotal)-min(freqtrTotal))/numxfit:max(freqtrTotal) ]';
    pfit{pp} = binomval_lims( b, xfit, 'probit', guessing, lapsing );
    plot(freqtrTotal, m.or, ['r*']);
    hold on
    plot( xfit, pfit{pp}, 'r', 'LineWidth', 2);
    set(gca, 'Ylim', [0 1]);
    [m.th50(pp)] = threshold_slope(pfit{pp},xfit,.5);
    [m.th20(pp), slope] = threshold_slope(pfit{pp},xfit,.20);
    [m.th80(pp), m.slope(pp)] = threshold_slope(pfit{pp},xfit,.80);
    
    plot(freqtrTotal, ones(1,length(freqtrTotal))*.5, '-.');
    set(gca, 'Ylim', [0 1], 'Fontsize', 15, 'Xlim', [freqtrTotal(1) freqtrTotal(end)]);
    plot([m.th20(pp) m.th20(pp)], [0 1], ':');
    plot([m.th80(pp) m.th80(pp)], [0 1], ':');
    x = xlabel('Pitch (Hz)', 'horizontalalignment','center');
    y = ylabel('Accuracy', 'horizontalalignment','center');
end

%% plot 2 example subjects and overall average
[r in] = sort(m.slope);
pptopl = [in(2) in(end-2)];
figure
for it = 1:length(pptopl)+1
    if it <= length(pptopl)
        subplot(1,3,it);
        pfit2 = pfit{pptopl(it)};
        pcor2 = pcor(pptopl(it),:)';
    elseif it == length(pptopl)+1
        subplot(1,3,3);
        pcor2 = pcor';
        pfit2 = cell2mat(pfit);
    end

    plotSpread(pcor2', 'xValues', freqtrTotal); hold on
    x = shadedErrorBar(xfit, mean(pfit2,2), std(pfit2,[],2)./sqrt(size(pfit2,2)));
    if it <= length(pptopl)
        plot(freqtrTotal, pcor2, ['r*']);
    else
        errorbar(freqtrTotal,mean(pcor2,2), std(pcor2,[],2)./sqrt(size(pcor2,2)), 'linestyle','none', 'linewidth', 2,'color', [1 0 0]);
    end
    
    x.patch.FaceAlpha =0.2;
    r = round(freqtrTotal([1 end]));
    set(gca, 'xtick', round(linspace(r(1), r(2),3)));
    
     set(gca, 'Ylim', [0 1]);
     plot(freqtrTotal, ones(1,length(freqtrTotal))*.5, 'k:');
     set(gca, 'Ylim', [0 1], 'Fontsize', 10, 'Xlim', [freqtrTotal(1) freqtrTotal(end)]);
     x = xlabel('Pitch (Hz)', 'horizontalalignment','center');
     y = ylabel('Accuracy', 'horizontalalignment','center');
end
set(gcf, 'position', [300 400 680 200]);




