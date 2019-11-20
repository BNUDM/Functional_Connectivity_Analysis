function result = ROIExtract(reference, roimask, transform, interpolation)
% ROIExtract Extract value from reference image within an ROI
%
% Usage: result = ROIExtract(reference, roimask [, transform] [, interpolation])
%
% See also: ROIData, ROICoord
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
if nargin < 4
    interpolation = false;
end

% Check and load reference image:
if ischar(reference)
    if ~exist(reference, 'file')
        error('%s: can not find file %s.', mfilename, reference);
    end
    featfolder = regexp(reference, ['.*\.feat', filesep], 'match', 'once');
    reference = ReadNii(reference);
    maskfile = [featfolder, 'mask.nii.gz'];
    if exist(maskfile, 'file')
        brainmask = ReadNii(maskfile);
        mask = brainmask.data ~= 0;
        reference.mask = mask;
    end
elseif ~isstruct(reference) || ~isfield(reference, 'sform')
    error('%s: Unrecognized input data %s.', mfilename, inputname(2));
end
if reference.info.dim ~= 3 && reference.info.dim ~= 4
    error('%s can only deal with 3D or 4D image.', mfilename);
end

% Check and load transform matrix:
if nargin >= 3
    identical = false;
    if ischar(transform)
        if ~exist(transform, 'file')
            error('%s: can not find file %s.', mfilename, transform);
        end
        transform = load(transform, '-ASCII');
    elseif isempty(transform)
        % identical, no transform used.
        identical = true;
    elseif ~ismatrix(transform) || size(transform,1) ~= 4 || size(transform,2) ~= 4
        error('%s: Unrecognized input matrix %s.', mfilename, inputname(3));
    end
else
    transform = [];
    identical = true;
end

if iscell(roimask)
    % for multiple ROIs:
    nrois = numel(roimask);
    result = cell(size(roimask));
    for ii = 1: nrois
        result{ii} = feval(mfilename, reference, roimask{ii}, transform, interpolation);
    end
    return;
elseif isstruct(roimask) && length(roimask) > 1
    % for multiple ROIs:
    nrois = numel(roimask);
%     result = repmat(struct('meants', [], 'ROI', '', 'TR', [], 't', [], 'timeseries', []), size(roimask));
    for ii = nrois: -1 : 1
        result(ii) = feval(mfilename, reference, roimask(ii), transform, interpolation);
    end
    return;
elseif ischar(roimask)
    if ~exist(roimask, 'file')
        error('%s: can not find file %s.', mfilename, roimask);
    end
    [~, roiname, ext] = fileparts(roimask);
    roimask = ReadNii(roimask);
    roimask.name = [roiname, ext];
elseif ~isstruct(roimask) || ~isfield(roimask, 'sform')
    error('%s: Unrecognized input data %s.', mfilename, inputname(1));
end

if roimask.info.dim ~= 3
    error('%s: ROI must be 3D image.', mfilename);
end

nvols = reference.info.nvolume;
if isfield(roimask, 'name')
    roiname = roimask.name;
else
    roiname = 'ROI';
end
if identical
    % check if the image are in the same space:
    if ~isequal(roimask.sform, reference.sform) || ...
            ~isequal(roimask.info.resolution, reference.info.resolution) || ...
            ~isequal(roimask.info.size, reference.info.size)
        error('%s: input reference and ROI were in different spaces.', mfilename);
    end
    % in the same space, so the data size should match:
    idx = roimask.data > 0.5;
    nvoxels = sum(idx(:));
    orig = zeros(nvols, nvoxels);
    for ii = 1: nvols
        tmp = reference.data(:,:,:, ii);
        orig(ii, :) = tmp(idx);
    end
else
    % perform the transform and obtain the coordinates:
    [maskidx, coords] = ROICoord(roimask, reference, transform);
    if interpolation
        % tri-linear interpolation to obtain the time series
        nvoxels = size(coords, 1);
        orig = zeros(nvols, nvoxels);
        for ii = 1: nvols
            tmp = interp3(reference.data(:,:,:, ii), coords(:,2), coords(:,1), coords(:,3));
            orig(ii, :) = tmp;
        end
    else
        % or use the rounded index
        idx1 = maskidx(:, 1) > 0 & maskidx(:, 1) <= roimask.info.size(1);
        idx2 = maskidx(:, 2) > 0 & maskidx(:, 2) <= roimask.info.size(2);
        idx3 = maskidx(:, 3) > 0 & maskidx(:, 3) <= roimask.info.size(3);
        idx  = maskidx(idx1 & idx2 & idx3, :);
        
        nvoxels = size(idx, 1);
        orig = zeros(nvols, nvoxels);
        for ii = 1: nvoxels
            orig(:, ii) = reference.data(idx(ii,1), idx(ii,2), idx(ii,3), :);
        end
    end
end % if identical transform matrix

if nvols > 2
    cor = corr(orig);
    % remove useless timeseries:
    emptyidx = isnan(diag(cor));
    if any(emptyidx)
        fprintf('Warning! %s contains %d empty voxels!', roiname, sum(emptyidx));
    end
    % orig(:, emptyidx) = [];
    % check the data quality:
    idx = ~triu(ones(nvoxels));
    cor = cor(idx);
    meancor = nanmean(cor);
    idx = cor < 0.3;
    lowerlevel = 100*mean(idx);
    lowercor = mean(cor(idx));
    if lowerlevel > 0
        fprintf('Warning! %s: mean correlation %.4f, %.2f%% lower than 0.3 with mean %.4f.\n', ...
            roiname, meancor, lowerlevel, lowercor);
    end
end

if ~isfield(reference, 'mask')
    % mask = sum(reference.data, 4) ~= 0;
    % because a correct input should have identical mask across t axis:
    mask = reference.data(:,:,:,1) ~= 0;
    mask = mask(:);
else
    mask = reference.mask(:);
end
[rr,cc,ss,vv] = size(reference.data);
tmp = reshape(reference.data, [rr*cc*ss, vv]);
temp = tmp(mask, :);
% because matlab's own mean is much slower, we'll calculate it manually
average = sum(temp, 1, 'double')/sum(mask);
result.average = average';

% meants is similar to fslmeants
meants = mean(orig, 2);
result.meants = meants;

result.ROI = roiname;
result.TR = reference.info.TR;
result.t  = reference.info.TR * (0.5: nvols)';
result.timeseries = orig;

end % function

% End Of File
