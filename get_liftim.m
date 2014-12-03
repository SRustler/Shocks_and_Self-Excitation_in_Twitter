%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% Notes:
% - This function finds the life time of the largest cluster (i.e. highest 
%   number of generations). 
% - Note that this time does NOT have to equal the maximal response time!
%
% Input:
% - cid: Cluster index
% - tevnts: time stamps of events
% 
% Output: 
% - h: Histogram of number of members per cluster
% - liftim: life time of the largest cluster
% 
function [h,liftim] = get_liftim(cid,tevnts)
    h = hist(cid,max(cid)); %Choose histogram such that each bin exactly corresponds to the number of members of that cluster.
    idxoflargestclu = find(h==max(h)); %Find indices of the largest clusters. Note that there can be several with the same maximal size!
    liftim = zeros(1,length(idxoflargestclu));
    for i = 1:length(idxoflargestclu)
        idx = find(cid==idxoflargestclu(i));
        liftim(i) = tevnts(idx(end)) - tevnts(idx(1));
    end
    liftim = max(liftim); %Out the of the largest cluster choose the more long-lived one.
end

%% This is how we calculated liftim earlier. 
% [h,inds] = sort(hist(cid,max(cid)),'descend'); %max(cid) is #clusters. hist checks how big each cluster is. The cluster index is arbitrary, thus sort. inds(1) will be index of largest cluster, inds
% h = hist(h,max(h)); %Just make histogram vector for plotting, overwriting the now obsolete h.
% idxoflargestclu = find(cid==inds(1)); %Vector of original indices that have the index of the largest cluster.
% liftim = tevnts(idxoflargestclu(end))-tevnts(idxoflargestclu(1)); %Note: liftim will be the largest response time, i.e. tau_max, which is relevant for the Laplace transformation: Int[...,{t,0,tau_max}] instead of Int[...,{t,0,infinity}]!
