function [result, hires] = ConvolveEV( ev, kernel, TR, volume )
% ConvolveEV Do convolve using a given kernel on the custom evfile.
%
% Usage: [result, hires] = ConvolveEV( ev, kernel, TR, volume )
%
% See also: conv, FslHRF
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
if ischar(ev)
    if ~exist(ev, 'file')
        error('Can not stat ev file %s.', ev);
    end
    evdata = load(ev, '-ASCII');
    evname = ev;
elseif ismatrix(ev)
    evdata = ev;
    evname = inputname(1);
else
    error('Unsupported input EV data %s.', input(1));
end

ncols = size(evdata, 2);
if ncols == 1
    % custom 1 column EV file
    if length(evdata) ~= volume
        error('No enough entries in EV %s.', evname);
    end
    onset = TR * (0: volume-1)';
    duration = TR * ones(volume, 1);
    weight = evdata;
elseif ncols == 3
    % custom 3 column EV file
    onset = evdata(:, 1);
    duration = evdata(:, 2);
    weight = evdata(:, 3);
else
    error('Unrecognized EV %s.', evname);
end

t = (0: kernel.dt: TR*volume)';
evweight = zeros(size(t));
to = 0;
for ii = 1: length(onset)
    tmp = onset(ii);
    idx = find(t>tmp, 1);
    if isempty(idx)
        warning('Out of range onset time %g. Skipped.', tmp);
        continue;
    end
    [delta1, adj] = min([t(idx)-tmp, tmp-t(idx-1)]);
    from = idx-(adj==2);
    if to >= from
        error('Onset time and duration is illegal in EV %s.', evname);
    end
    tmp = onset(ii) + duration(ii);
    idx = find(t>tmp, 1);
    if isempty(idx)
        warning('Out of range event time %g. Skipped.', tmp);
        continue;
    end
    [delta2, adj] = min([t(idx)-tmp, tmp-t(idx-1)]);
    to = idx-(adj==2);
    if delta1 >= kernel.dt || delta2 >= kernel.dt
        warning('Onset time or duration might be out of range in %s.', evname);
    end
    evweight(from: to) = weight(ii);
end

% FSL did a `remmean` before convolution, and it seems useless, 
% except that they have to append some `negpts` to workaround it.
% evweight = evweight - mean(evweight);

kern = kernel.kernel ./ sum(kernel.kernel);
evmodel = conv(evweight, kern, 'full');
hires = evmodel(1: length(t));
% NOTE that FSL's `do_convolve` is buggy, but because of the same workaround
% I've mentioned above, this bug does not cause too much trouble.

% down-sampling: FSL use a "MIDDLE" (t+0.5) for resample down,
% Considering about slice timing correction, we use the real middle here:
result = hires( 1+round((0.5: volume)*TR/kernel.dt) );

% Finally, FSL did a `peakAndFilter`, which is not performed here.
end % function

% End Of File
