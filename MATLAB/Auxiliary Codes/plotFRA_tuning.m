function median_post = plotFRA_tuning(FRA, myKsDir, clusters, n_bins)
%
% Function that prints a plot with the histograms of the spikes prior and
% post trigger per frequency in the FRA of the chosen clusters (log-scale).
% All the information is contained in the structure FRA that is yielded 
% from the function clusterInFRA.
% 
% Inputs:
%     FRA: struct, contains all the fields pertaining to FRA sequence.
%     myKsDir: str, where the KiloSort files are located.
%     clusters: clusters we want to plot.
%     n_bins: number of bins in which the frequencies are divided.
%   
% Outputs:
%     median_post: double, vector with the median response frequency of
%     each cluster.

% Load variables.
num_cl = length(clusters);
freqs = FRA.freqTrig;

pre_trig = FRA.pre_trig;
post_trig = FRA.post_trig;

% Figure out the disposition of the graph in rows and cols.
rows = floor(sqrt(num_cl));
cols = floor(num_cl/rows) + sum(mod(num_cl,rows)~=0);
    
median_post = [];

% Print the graphs.
for k=1:num_cl
    
    % Not making subplots when we have just 1 cluster.
    if num_cl >1
        subplot(rows,cols,k)
    end
    
    % Sort the frequencies in order so we can do an histogram with them.
    [a,b]=sort(freqs);
    spikes_post =post_trig(b,clusters(k)+1);

    num_spikes_post = [];
    for i=1:length(a)
        if spikes_post(i)~=0
            for j=1:spikes_post(i)
                num_spikes_post = [num_spikes_post, a(i)]; %#ok<*AGROW>
            end
        end
    end
    
    [a,b]=sort(freqs);
    spikes_pre =pre_trig(b,clusters(k)+1);
    
    % Calculating the number of spikes that exist and their associated
    % frequencies.
    num_spikes_pre = [];
    for i=1:length(a)
        if spikes_pre(i)~=0
            for j=1:spikes_pre(i)
                num_spikes_pre = [num_spikes_pre, a(i)];
            end
        end
    end
    
    median_post(k) = 0;
    
    % Calculating the medians.
    if ~isempty(num_spikes_post)
        median_post(k) = median(num_spikes_post);
        median_to_show = num2str(round(median_post(k)/1000,1));
    else 
        median_to_show = "No median";
    end
   
    % Creating bins and counts in logarithmic scale.
    edges = 10.^(3:2.5/(n_bins):(5.5-0.001));
    [N_pre, edges_pre] = histcounts(num_spikes_pre,edges);
    [N_post, edges_post] = histcounts(num_spikes_post,edges);
    
    % Print the histograms.
    hold on
    
    h(1) = histogram('BinEdges',edges_post,'BinCounts',N_post,'DisplayName','Post trig');
    h(2) = histogram('BinEdges',edges_pre,'BinCounts',N_pre,'DisplayName','Pre trig');
    h(3) = xline(median_post(k), '--r','DisplayName','Median');
    
    % Set y limits and a logarithmic scale in x.
    set(gca, 'Xscale',  'log');
    ylim([0, (max([N_pre N_post])*1.2)])
    
    % Add labeling and legends.
    [amp, good] = findAmpAndSorting(myKsDir, clusters(k));
    title(clusters(k) + ", " + good + ", " + amp + " mV, "+ median_to_show+ " kHz")
    
    % Only show the legend at the last cluster (no need at all of them).
    if k == num_cl
        legend(h(1:3))
    end
    
    % Axes labelling.
    ylabel('# of spikes')
    xlabel('Frequency (Hz)')
end