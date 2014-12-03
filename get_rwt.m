%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% Notes: 
% - This function takes the information on the cluster of each event to 
%   calculate the renormalized waiting times (rwt) and much more. See Input
%   and Output below.
% 
% Input:
% - CID: Cluster index of each event
% - t: Sorted time stamps of all events
% - gens: Associated generation of time stamps
% - threshold: Lower threshold for cluster size to be considered for
%   rwt-calculation. Clusters smaller than threshold will be discarded.
%   This corresponds to the synthetic case of examining only the "Top-
%   tweets", which are those that are re-tweeted most (among other 
%   criteria). Events of clusters that fall below the threshold will be 
%   marked by -1. Default value = 0.
% - recents: Number of most recent responses to be considered for
%   rwt-calculation. Responses older than recents will be marked by -1.
%   Default value = nan.
% 
% Output:
% - rwt: Renormalized waiting times, i.e. the waiting time of each event to 
%   its respective cluster center, i.e. immigrant. Events not fulfilling
%   criteria specified by threshold and recents, will be marked by -1.
% - K: Matrix with each row corresponding to the independent clusters and
%   the n-th column corresponding to the longest time it took for the k-th
%   generation to occur since immigration.
% - t0: Time stamps of cluster centers (Recall: They follow a homogeneous 
%   Poisson process)
%
function [rwt,K,t0] = get_rwt(CID,t,gens,threshold,recents)
    if nargin<5
        recents = nan;
    end

    if nargin<4
        threshold = 0;
    end

    [cid,inds] = sort(CID); %Default: ascending such that cid(i) <= cid(i+1)
    t = t(inds); %Sort t in the same way CID was sorted.
    gens = gens(inds);
    maxgens = max(gens);
    maxcid = max(CID); %Index of latest cluster to occur.

    %Find cluster borders
    % Logic:
    % i:     1     2     3     4     5     6     7     8     9    10   11 ... length(t)
    % cid:   1     1     1     2     2     2     2     3     3     3    4 ... length(t)
    % first: 1                 4                       8               11 ... maxcid
    % last:              3                       7                10      ... maxcid
    % => first(N) and last(N) denote index (sorted via cid) of the first 
    %    and last time stamp of the N-th cluster.
    
    ccount = 0; %Cluster counter
    last = zeros(1,maxcid); %For each cluster we want to know the index of its LAST time stamp.
    h = find(cid==maxcid); %Store indices of events of the latest cluster to occur. 
    last(maxcid) = h(end); %Manually assign value because otherwise it would be left at 0. We know that the LAST event of the last cluster (to occur) has index h(end).
    first = zeros(1,maxcid); %For each cluster we want to know the index of its FIRST time stamp.
    first(1) = 1; %Again, manually assign value. We know that the FIRST event of the first cluster (to occur) was at t(1).
    for i=1:length(cid)-1 %Loop through all events sorted via cid
        if cid(i) < cid(i+1) %Cluster border, i.e. new cluster, found.
            ccount = ccount + 1; %Increase cluster count

            last(ccount) = i; %Set "index of previously found cluster's LAST time stamp" to current time stamp.
            first(ccount+1) = i+1; %Set "index of newly found cluster's FIRST time stamp" to next time stamp.
        end
    end

    %Find waiting times to immigrants
    t0 = zeros(1,maxcid);
    rwt = zeros(1,length(t)); %Renormalized waiting times
    K = zeros(maxcid,maxgens);
    for i=1:maxcid  %Loop through all clusters
        t0(i) = min(t(first(i):last(i))); %E.g. min(t[4 5 6 7]): the earliest time stamp of the second cluster is our second t0, i.e. t0(2).
        temp = t(first(i):last(i)) - t0(i); %Calculte waiting times of the i-th cluster renormalized to its cluster center having occured at t0(i).
        rwt(first(i):last(i)) = temp;
        
        if length(temp)>recents %If the cluster size is higher than recents...
            rwt(first(i):last(i)-recents) = -1; %...assign -1 to earliest responses s.t. only recents most recent responses remain. This is the synthetic case of being only able to get the last 100 RTs of a MT.
        end
        if length(temp)<threshold %If the cluster size falls below a threshold...
            rwt(first(i):last(i)) = -1; %...assign -1 to be filtered out. rwt(rwt>-1) will be ALL events fulfilling the threshold. This is needed for the calculation of n'. This is the synthetic case examining the "TOP-Tweets" only.
        end
        
        for k = 1:max(gens(first(i):last(i)))
            % If k empty KK will not be affected.
            K(i,k) = max(temp(gens(first(i):last(i)) == k)); %Check when gens of current cluster equal current k, then take only the k-event with the highest rwt (s.t. tmax more likely to be in it)
        end
        
    end
    t0 = sort(t0); %Confirmed: t0 is exactly texo!
    [~,inds] = sort(t);
    rwt = rwt(inds); %Sort rwt in the same way t was sorted such that the j-th entry of rwt corresponds to the j-th entry of t.
end
