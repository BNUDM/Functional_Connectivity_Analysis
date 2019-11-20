function [result, data] = DeconvGLM( config, data )
% DeconvGLM deconvolve time series into trial-by-trial beta value using GLM
%
% Usage: [result, data] = DeconvGLM(config [, data])
%           where config is an INI file, see deconv.cfg for details.
%
% See also: MeanTS
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
ticstart = tic;
if nargin == 0
    [fname, fpath] = uigetfile({'*.ini;*.cfg;*.conf', 'Configure file'}, ...
        'Please select the timeseries config file');
    config = fullfile(fpath, fname);
    fprintf('>> %s(''%s'')\n', mfilename, config);
end

if ~exist(config, 'file')
    error('Can not stat configure file %s.', config);
end

%%
% configure file is read by ROIData
fprintf('>>> Stage 0, reading data...\n');
if nargin == 2
    % use the input argument instead of reading cfg file
    [~, cfg, option, info] = ROIData(config, 'dryrun', true);
    % check compatibility
    [dm, dn] = size(data);
    ntasks = length(info.task);
    nrois  = length(cfg.stdroi);
    if dm ~= ntasks || dn ~= nrois
        warning('Input data is not compatible with config file!');
    end
else
    [data, cfg, option, info] = ROIData(config);
end

%%
% Read EV data:
[event, option] = EVData(cfg, option, info);

% convolution kernel: canonical HRF
hrf = FslHRF(0: option.interval: 36);
%%
% for each task, split into segments and transform into a single number
ntasks = length(info.task);
nrois  = length(cfg.stdroi);
nsubjects = length(info.subjects);

% initialize result:
for ii = 1: size(event, 2)
    tmp.(event(1, ii).condition).(event(1, ii).stage) = [];
end
temp = tmp;
result = repmat(tmp, nsubjects, nrois);

fprintf('>>> Stage 1, fitting and concatenating...\n');
for ii = 1: ntasks
%     fprintf('.');
    fprintf('[%3d%%] Processing on %s...\n', round(ii/ntasks*100), info.task{ii});

    nvolumn = length(data(ii,1).meants);
    TR = data(ii,1).TR;
    
    designmat = [];
    evidx = temp;
    for kk = 1: size(event, 2)
        ev = event(ii, kk);
        ntrial = size(ev.mat,1);
        evmat = zeros(nvolumn, ntrial);
        for jj = 1: ntrial
            evmat(:,jj) = ConvolveEV(ev.mat(jj,:), hrf, TR, nvolumn);
        end
        evidx.(ev.condition).(ev.stage) = size(designmat,2)+(1: ntrial);
        designmat = cat(2, designmat, evmat);
    end % for all events
    
    designmat = bsxfun(@minus, designmat, mean(designmat));
    confound = [data(ii,1).average, data(ii,1).mc];
    confound = bsxfun(@minus, confound, mean(confound));
    
    % matrix version of LM should be much faster:
    realconf = [ones(size(confound,1), 1), confound];
    iconf    = pinv(realconf);
    idmat    = pinv(designmat);
    
    for jj = 1: nrois
        if option.voxelwise
            tmp = double(data(ii,jj).timeseries);
        else
            tmp = double(data(ii,jj).meants);
        end
        bconf = iconf * tmp;
        dataf = realconf * bconf;
        resid = tmp - dataf;
        beta = idmat * resid;
        
        for kk = 1: size(event, 2)
            ev = event(ii, kk);
            result(info.index(ii), jj).(ev.condition).(ev.stage) = ...
                cat(1, result(info.index(ii), jj).(ev.condition).(ev.stage), ...
                beta(evidx.(ev.condition).(ev.stage), :));
        end % for each event
    end % for each roi
end % for each task

fprintf('All done in %gs.\n\n', toc(ticstart));

end % function

% End Of File
