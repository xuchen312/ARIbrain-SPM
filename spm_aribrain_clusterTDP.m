function [tdpclus] = spm_aribrain_clusterTDP(allp,sortp,ixp,varargin)
%
% Compute the ARI-based TDP lower bounds for any clusters
%
% =========================================================================
% FORMAT: [tdpclus] = spm_aribrain_clusterTDP(allp,sortp,ixp,['alpha', ...
%                     alpha,'simes',simes])
% -------------------------------------------------------------------------
% Inputs:
%  -    allp: a vector of all p-values
%  -   sortp: a vector of all sorted p-values in ascending order
%  -     ixp: a cell array of 1-based index subsets for all clusters
% Inputs (optional):
%  -   alpha: a significance level (default: 0.05)
%  -   simes: a logical variable used to indicate Simes' (true, default) or 
%             Hommel's (false) robust test is conducted
%
% Outputs:
%  - tdpclus: a vector of TDP lower bounds for all input clusters
% =========================================================================
%

% ------- (1) Initialization & preparation

% Check inputs & outputs
if mod(nargin,2) ~= 1 
    error('spm_aribrain_clusterTDP: there must be an odd number of inputs'); 
end
if nargin < 3; error('spm_aribrain_clusterTDP: not enough inputs'); end
if nargin > 7
    error('spm_aribrain_clusterTDP: too many input arguments'); 
end
if nargout ~= 1 
    error('spm_aribrain_clusterTDP: there must be one output'); 
end

% Specify default optional input values
alpha = 0.05;
simes = true;
% Load optional input arguments
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'alpha'
            alpha = varargin{2};
        case 'simes'
            simes = varargin{2};
        case 'conn'
            warning('spm_aribrain_clusterTDP: ''conn'' not used');
        case 'tdp'
            warning('spm_aribrain_clusterTDP: ''tdp'' not used');
        otherwise
            error(['spm_aribrain_clusterTDP: unexpected argument: ' varargin{1}])
    end
    varargin(1:2) = [];
end

% ------- (2) Estimate the TDP

NUM = length(ixp);  % cluster number

% loop through clusters
tdpclus = zeros(NUM,1);
for i=1:NUM
    tdpclus(i) = double(spm_aribrain_cluster(length(allp),length(ixp{i}), ...
        alpha,simes,int32(ixp{i}),allp,sortp))/length(ixp{i});
end

end
