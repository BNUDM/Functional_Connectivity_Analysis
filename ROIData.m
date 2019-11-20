function [data, cfg, option, info] = ROIData(config, varargin)
% ROIData Extract data within ROIs from multiple feat folders
%
% Usage: [data, cfg, option, info] = ROIData(config)
%
% See also: ROIExtract
%
% This function is released as part of the `FSLboost` package.
%                                   Powered by Matlab && FSL.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2014-2016 LiTuX @BNU, all rights reserved.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if ~exist(config, 'file')
    error('Can not stat configure file %s.', config);
end

%%
% config is a kind of ini file, read it and parse:
cfg = ReadINI(config);
rpath = GetCanonicalPath(config);

% Global section: define input et al.
tmp = struct2cell(cfg.global);
option = cell2struct(tmp(2,:)', tmp(1,:)');

if isfield(option, 'interpolation') && strcmpi(option.interpolation, 'true')
    option.interpolation = true;
else
    option.interpolation = false;
end

if ~isfield(option, 'input')
    option.input = 'filtered_func_data.nii.gz';
end

%%
dryrun = false;
% overwritten option using input arguments
for ii = 1: 2: (nargin-1)
    if isfield(option, varargin{ii})
        option.(varargin{ii}) = varargin{ii+1};
    elseif strcmpi(varargin{ii}, 'dryrun')
        dryrun = varargin{ii+1};
    else
        warning('Unknown option %s.', varargin{ii});
    end
end

%%
% read the input data information:
if IsAbsolutePath(option.subjectinfo)
    sfile = option.subjectinfo;
else
    sfile = fullfile(rpath, option.subjectinfo);
end
if IsAbsolutePath(option.taskinfo)
    tfile = option.taskinfo;
else
    tfile = fullfile(rpath, option.taskinfo);
end
info = DataInfo(sfile, tfile);

%%
if dryrun
    data = [];
    return;
end

%%
% stdroi section: pre load the roi mask
nrois = length(cfg.stdroi);
for ii = nrois: -1: 1
    tmp = cfg.stdroi(ii).value;
    if IsAbsolutePath(tmp)
        roi = ReadNii(tmp);
    else
        % relative to the config file
        roi = ReadNii(fullfile(rpath, tmp));
    end
    roi.name = cfg.stdroi(ii).key;
    roimask(ii) = roi;
end

%%
% for each task, extract the mean timeseries:
ntasks = length(info.task);

for ii = ntasks: -1: 1
    tmp = [option.suffix, '.feat'];
    
    epifile = [info.task{ii}, tmp, filesep, option.input];
    if ~exist(epifile, 'file')
        error('%s: can not stat epi image for %s.', mfilename, info.task{ii});
    end
    transform = [info.task{ii}, tmp, filesep, 'reg', filesep, 'standard2example_func.mat'];
    if ~exist(transform, 'file')
        error('%s: can not stat transform for %s.', mfilename, info.task{ii});
    end
    
%     mcfile = [info.task{ii}, tmp, filesep, 'mc', filesep, 'prefiltered_func_data_mcf_final.par'];
    mcfile = [info.task{ii}, tmp, filesep, 'mc', filesep, 'prefiltered_func_data_mcf.par'];
    mcdata = load(mcfile, '-ASCII');
    
    fprintf('[%3d%%] Reading %s.\n', round(100-ii/ntasks*100), epifile);
    ts = ROIExtract(epifile, roimask, transform, option.interpolation);
    
    [ts.mc] = deal(mcdata);
    
    data(ii,:) = ts;
end

end % function

