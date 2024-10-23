function [nclus,Lclus,tdpclus] = spm_aribrain_tdpCluster(allp,sortp,ordp,rankp,mask,varargin)
%
% Compute maximal clusters using adaptive thresholding algorithm
%
% =========================================================================
% SYNTAX: [nclus,Lclus,tdnclus] = spm_aribrain_tdpCluster(allp,sortp,ordp, ...
%                                 rankp,mask,['alpha',alpha,'conn',conn, ...
%                                 'simes',simes,'tdpth',tdpth])
% -------------------------------------------------------------------------
% Inputs:
%  -    allp: a vector of all valid p-values
%  -   sortp: a vector of all sorted p-values in ascending order
%  -    ordp: a vector of sort indices, which rearranges allp to sortp
%  -   rankp: a vector of sorting ranks
%  -    mask: a 3D array, where zero indicates unwanted voxels
% Inputs (optional):
%  -   alpha: a significance level (default: 0.05)
%  -    conn: a connectivity criterion (default: 18)
%  -   simes: a logical variable used to indicate Simes' (true, default) or 
%             Hommel's (false) robust test is conducted
%  -   tdpth: a chosen TDP threshold (default: 0.7)
%
% Outputs:
%  -   nclus: a vector of generated maximal cluster sizes
%  -   Lclus: a 3D array of clusters, labelled with 1:#{clusters}
%  - tdpclus: a vector of TDP lower bounds for all maximal clusters
% =========================================================================
%

% ------- (1) Preparations

% check inputs & outputs
if mod(nargin,2) ~= 1 
    error('spm_aribrain_tdpCluster: there must be an odd number of inputs'); 
end
if nargin < 5; error('spm_aribrain_tdpCluster: not enough inputs'); end
if nargin > 13; error('spm_aribrain_tdpCluster: too many inputs'); end
if nargout ~= 3; error('spm_aribrain_tdpCluster: there must be 3 outputs'); end

% specify default optional input values
alpha = 0.05;
conn  = 18;
simes = true;
tdpth = 0.7;

% load optional input arguments
while ~isempty(varargin)
    switch lower(varargin{1})
        case 'alpha'
            alpha = varargin{2};
        case 'conn'
            conn = varargin{2};
        case 'simes'
            simes = varargin{2};
        case 'tdpth'
            tdpth = varargin{2};
        otherwise
            error(['spm_aribrain_tdpCluster: unexpected argument: ' varargin{1}])
    end
    varargin(1:2) = [];
end

% ------- (2) Form clusters

indexp = find(mask~=0);   % in-mask voxel indices in 3D space
dims   = size(mask);      % image dimensions
m      = length(allp);    % number of all p-values
maskI  = zeros(dims);
maskI(indexp) = 1:m;      % 3D mask of original orders (1:m)

% create maximal clusters based on tdpth
[nclus,ixclus,tdpclus] = spm_aribrain_tdp(m,conn,alpha,tdpth,simes,int32(dims), ...
    int32(maskI(:)),int32(indexp),int32(ordp),int32(rankp),allp,sortp);

% sort clusters in descending order by size
NUM = length(nclus);      % cluster number
[nclus,ordc] = sort(nclus,'descend');
rankc = zeros(NUM,1);     % get sorting ranks
rankc(ordc) = 1:NUM;
tdpclus = tdpclus(ordc);  % sort TDP bounds

% compute Lclus
Lclus = zeros(dims);
start = 1;
for i=1:NUM
    last = sum(nclus(1:i));
    Lclus(indexp(ixclus(start:last)+1)) = rankc(i);
    start = last+1;
end

end
