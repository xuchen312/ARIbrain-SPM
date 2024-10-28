function varargout = spm_aribrain(varargin)
%
% Run TDP inference using ARI
%
% =========================================================================
% FORMAT:                          spm_aribrain(['xSPM',xSPM,'alpha',alpha,'file',file,'simes',simes,'conn',conn,'tdpth',tdpth])
%         [hReg,xSPM,SPM,TabDat] = spm_aribrain(['xSPM',xSPM,'alpha',alpha,'file',file,'simes',simes,'conn',conn,'tdpth',tdpth])
% -------------------------------------------------------------------------
% Inputs (optional; if empty or not specified, the default is used):
%  -    xSPM: an input structure containing SPM, distribution & filtering
%             details (see spm_getSPM.m for details; default: compute xSPM
%             interactively)
%  -   alpha: a significance level (default: 0.05)
%  -    file: a character array specifying the file name for an output CSV
%             file (default: output table is not saved)
%  -   simes: a logical variable used to indicate Simes' (true, default) or
%             Hommel's (false) robust test is conducted
%  -    conn: a connectivity criterion; 6 (face), 18 (edge, default) or 26 
%             (vertex) (only needed for the tdpCluster inference)
%  -   tdpth: a chosen TDP threshold (default: the clusterTDP inference 
%             will be performed; only needed for the tdpCluster inference)
%
% Outputs (optional, for interactive exploration):
%  -    hReg: handle of MIP XYZ registry object
%             (see spm_XYZreg.m for details)
%  -    xSPM: an evaluated/thresholded structure containing SPM,
%             distribution & filtering details
%             (see spm_getSPM.m for details)
%  -     SPM: an SPM structure
%             (see spm_getSPM.m for details)
%  -  TabDat: a structure containing table data with fields
%             (see spm_clusterTDP_list.m for details)
% =========================================================================
%

%-Set default modality
%----------------------------------------------------------------------
try
    modality = spm_get_defaults('modality');
    spm('Defaults',modality);
catch
    spm('Defaults','FMRI');
end

%-Check input arguments
%----------------------------------------------------------------------
if mod(nargin,2) ~= 0
    error('spm_aribrain: there must be an even number of input arguments');
end
if nargin > 12; error('spm_aribrain: too many input arguments'); end
% load inputs
while ~isempty(varargin)
    switch varargin{1}
        case 'xSPM'
            if ~isempty(varargin{2}); xSPM = varargin{2}; end
        case 'alpha'
            if ~isempty(varargin{2})
                alpha = varargin{2};
                if ~isnumeric(alpha) || alpha <=0 || alpha >=1
                    error('spm_aribrain: ''alpha'' must be numeric & within (0,1)');
                end
            end
        case 'file'
            if ~isempty(varargin{2})
                file = varargin{2};
                if ~ischar(file)
                    error('spm_aribrain: ''file'' must be a character array');
                else
                    [fpath,fname,fext] = fileparts(file);
                    if isempty(fext)
                        file = strcat(file,'.csv');
                    elseif ~strcmpi(fext,'.csv')
                        error('spm_aribrain: unexpected output file extension: %s',fext);
                    end
                    if isempty(fname)
                        file = fullfile(fpath,strcat('aribrain',fext));
                        warning('spm_aribrain: found empty file name and use ''aribrain.csv''');
                    end
                end
            end
        case 'simes'
            if ~isempty(varargin{2})
                simes = varargin{2};
                if ~islogical(simes)
                    error('spm_aribrain: ''simes'' must be logical');
                end
            end
        case 'conn'
            if ~isempty(varargin{2})
                conn = varargin{2};
                if conn ~= 6 && conn ~= 18 && conn ~= 26
                    error('spm_aribrain: ''conn'' must be 6 (face), 18 (edge) or 26 (vertex)');
                end
            end
        case 'tdpth'
            if ~isempty(varargin{2})
                tdpth = varargin{2};
                if ~isnumeric(tdpth) || tdpth < 0 || tdpth > 1
                    error('spm_aribrain: ''tdpth'' must be numeric & within [0,1]');
                end
            end
        otherwise
            error(['spm_aribrain: unexpected argument: ' varargin{1}])
    end
    varargin(1:2) = [];
end

%-Query SPM and setup GUI
%----------------------------------------------------------------------
if exist('xSPM','var')
    [hReg,xSPM,SPM] = spm_aribrain_ui('Setup',xSPM);
else
    [hReg,xSPM,SPM] = spm_aribrain_ui('Setup');
end

%-Compute TDP estimation summary table "TabDat"
%----------------------------------------------------------------------
if exist('alpha','var'); xSPM.alpha = alpha; end
if exist('simes','var'); xSPM.simes = simes; end
if exist('conn', 'var'); xSPM.conn  = conn;  end
if exist('tdpth','var'); xSPM.tdpth = tdpth; end
TabDat = spm_aribrain_list('List',xSPM,hReg);

%-Display table "TabDat"
%----------------------------------------------------------------------
spm_aribrain_list('TxtList',TabDat);

%-Write TDP estimation result table to a csv file
%----------------------------------------------------------------------
if exist('file','var')
    spm_aribrain_list('CSVList',TabDat,file);
end

%-Return outputs for interactive exploration in control panel
%----------------------------------------------------------------------
if nargout > 0
    varargout = {hReg,xSPM,SPM,TabDat};
end

return
