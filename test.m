% ------- test spm_aribrain_clusterTDP.m
Tmap = spm_read_vols(spm_vol('/Users/xuchen/Documents/work/LUMC/spm_ARI/R/livio/GOnoGO/z.nii'));
Pmap = spm_Ncdf(-Tmap);  % one-sided p-values
Pmap = spm_read_vols(spm_vol('/Users/xuchen/Documents/work/LUMC/spm_ARI/R/livio/GOnoGO/p.nii'));
mask = Tmap~=0;
allp = Pmap(mask~=0);

Tmap = spm_read_vols(spm_vol('/Users/xuchen/Documents/clusterTDP/lower_bound/SPM_extension/auditory/spmT_0001.nii'));
Pmap = spm_Tcdf(-abs(Tmap),xSPM.df(2))*2;  % two-sided p-values
mask = Tmap~=0;
allp = Pmap(mask~=0);

[sortp,~] = sort(allp);
Tmap(Tmap<3.1) = 0;
[L,NUM] = spm_bwlabel(Tmap,18);
iL = L(mask~=0);

tic;
ixp = cell(NUM,1);
for i=1:NUM
    ixp{i} = find(iL(:)==i);
end
varargins = {'alpha',0.01,'simes',false};
tdpclus = spm_aribrain_clusterTDP(allp,sortp,ixp,varargins{:})
T1=toc

tic;
tdpclus = zeros(NUM,1);
for i=1:NUM
    ixp = find(iL(:)==i);
    tdpclus(i) = spm_aribrain_clusterTDP(allp,sortp,{ixp},'simes',false);
end
tdpclus
T2=toc

[hReg,xSPM,SPM,TabDat] = spm_aribrain;
[hReg,xSPM,SPM,TabDat] = spm_aribrain('conn',26);
[hReg,xSPM,SPM,TabDat] = spm_aribrain('simes',false);
[hReg,xSPM,SPM,TabDat] = spm_aribrain('conn',[],'file','aribrain');
[hReg,xSPM,SPM,TabDat] = spm_aribrain('alpha',0.01,'simes',false);

% ------- test spm_aribrain_tdpCluster.m
Tmap = spm_read_vols(spm_vol('/Users/xuchen/Documents/work/LUMC/spm_ARI/R/livio/GOnoGO/z.nii'));
Pmap = spm_Ncdf(-Tmap);
Pmap = spm_read_vols(spm_vol('/Users/xuchen/Documents/work/LUMC/spm_ARI/R/livio/GOnoGO/p.nii'));
mask = Tmap~=0;
allp = Pmap(mask~=0);
m = length(allp);

Tmap = spm_read_vols(spm_vol('/Users/xuchen/Documents/clusterTDP/lower_bound/SPM_extension/auditory/spmT_0001.nii'));
Pmap = spm_Tcdf(-abs(Tmap),xSPM.df(2))*2;  % two-sided p-values
mask = Tmap~=0;
allp = Pmap(mask~=0);
m = length(allp);

mask = mask.*(~ismissing(Tmap));

[sortp,ordp] = sort(allp);
rankp        = zeros(m,1);
rankp(ordp)  = 1:m;

tic;
[nclus,Lclus,tdpclus] = spm_aribrain_tdpCluster(allp,sortp,ordp,rankp,mask);
T3=toc

[hReg,xSPM,SPM,TabDat] = spm_aribrain('tdpth',0.7)




% ------- some relevant code

% extract data from SPM
dims = SPM.xVol.DIM;   % image dimensions
STAT = SPM.xCon.STAT;  % statistic indicator character ('T', 'F' or 'P')
df   = SPM.xX.erdf;    % effective residual degrees of freedom
%df   = [1 SPM.xX.erdf];  % [df{interest} df{error}]

% cdtth: if within [0,1], assume it to be a p-value
cdtth = spm_input(['Cluster defining threshold {',STAT,' or p-value}'],'+0','r',0.001,1);
if (~isempty(cdtth) && cdtth >= 0 && cdtth <= 1)
    cdtth = spm_u(cdtth,[1 df],STAT);
end
% tdpth: within [0,1]
tdpth = spm_input('TDP threshold','+0','r',0.7,1,[1,0]);
% output file name: output must be a csv file
file = spm_input('Output CSV filename','+0','s','aribrain.csv');

% convert statistics to p-values
switch STAT
    case 'F'
        Pmap = 1-spm_Fcdf(Tmap,df);
    case 'T'
        Pmap = spm_Tcdf(-abs(Tmap),df)*2;  % two-sided p-values
        %Pmap = 1-spm_Tcdf(Tmap,df);   % one-sided p-values
    case 'X'
        Pmap = spm_Xcdf(Tmap,df,'upper');
    case 'Z'
        Pmap = spm_Ncdf(-abs(Tmap))*2;     % two-sided p-values
        %Pmap = 1-spm_Ncdf(Tmap);      % one-sided p-values
end

% extract SPM image from SPM
if (isfield(SPM,'xCon') && isfield(SPM.xCon,'Vspm') && ~isempty(SPM.xCon.Vspm))
    Tmap = spm_read_vols(spm_vol(SPM.xCon.Vspm));  % there could be more than one contrasts
else
    error('spm_aribrain: cannot find contrast-related SPM images');
end

% check mask
if (isfield(SPM,'VM') && ~isempty(SPM.VM))
    mask = spm_read_vols(spm_vol(SPM.VM));
    if ~isequal(size(Tmap), size(mask))
        error('spm_aribrain: SPM image & mask must have the same dimensions');
    else
        mask = mask.*(~ismissing(Tmap));  % discard voxels with missing values
    end
else
    warning('spm_aribrain: mask is created by removing missing values');
    mask = ~ismissing(Tmap);
end
if sum(mask~=0)==0
    error('spm_aribrain: mask contains zero valid voxels');
end

% cluster-forming procedure
Tmap(Tmap<cdtth) = 0;
[L,NUM] = spm_bwlabel(Tmap,18);
iL = L(mask~=0);

ixp = cell(NUM,1);
for i=1:NUM
    ixp{i} = find(iL(:)==i);
end


