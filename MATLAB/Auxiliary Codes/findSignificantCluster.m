function [cl_final, pre_trig, post_trig, p] = findSignificantCluster(spTimes, spSites, cl_chosen, triggersInterest, sec_before_trigger)
%
% Finding good clusters as the ones with significant difference between
% before and after the trigger.
% 
% Inputs:
%     spTimes: vector with doubles, times at which the spikes occurred.
%     spSites: vector with doubles, clusters to which the spikes pertain.
%     cl_chosen: num vector, all clusters you want to analyze.
%     triggersInterest: triggers appearing in the SOI.
%     sec_before_trigger: num, how many seconds before the triggers we want
%     to see to perform the t-test.
%
% Outputs:
%     cl_final: int, vector containing the clusters that have significant
%     differences.
%     pre_trig: matrix, number of triggers x # of clusters, number of
%     spikes of that cluster just prior to that trigger.
%     post_trig: matrix, number of triggers x # of clusters, number of
%     spikes of that cluster just after that trigger.
%     p: double, vector containing the p-values of the t-test compairing 
%     the difference prior and after the trigger.

% Initializing variables.
pre_trig = zeros(length(triggersInterest), length(cl_chosen));
post_trig = zeros(length(triggersInterest), length(cl_chosen));
h = zeros(length(cl_chosen),1);
p = zeros(length(cl_chosen),1);

for i = 1: length(cl_chosen)
    % Find which spikes correspond to that cluster.
    trial_sp = spTimes(spSites==cl_chosen(i))';
    
    % Finding the number of spikes that are contained before and after 
    % each trigger (integrating over t=trigger to t=trigger +-
    % sec_before_trigger).
    for k = 1:length(triggersInterest)
        pre_trig(k,i) = sum(trial_sp<triggersInterest(k) & trial_sp>(triggersInterest(k)-sec_before_trigger));        
        post_trig(k,i) = sum(trial_sp>triggersInterest(k) & trial_sp < (triggersInterest(k)+sec_before_trigger));
    end
    
    % Performing the t-test (whether there is a significant difference on
    % the amount of spikes before and after the triggers.
    [h(i), p(i)] = ttest2(pre_trig(:,i),post_trig(:,i));
end

% Storing all the clusters that reject the null hypothesis.
cl_final = cl_chosen(h==1);