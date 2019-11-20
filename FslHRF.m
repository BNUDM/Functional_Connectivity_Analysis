function kernel = FslHRF( t, delay1, sigma1, delay2, sigma2, ratio )
% FslHRF: return the double gamma HRF used by FSL & SPM et al.
%
% Usage: kernel = FslHRF( t, delay1, sigma1, delay2, sigma2, ratio )
%
% all arguments are optional with default value (same as in FSL):
%       t = 0: 0.05: 36;    ratio = 6;
%       delay1 = 6;         sigma1 = sqrt(6);
%       delay2 = 16;        sigma2 = sqrt(16);
% see `fsl/src/feat5/feat_model.cc` for details.
%
% See also: gampdf
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

if nargin < 1
    t = [];
end
if nargin < 2
    delay1 = 6;
end
if nargin < 3
    sigma1 = [];
end
if nargin < 4
    delay2 = 16;
end
if nargin < 5
    sigma2 = [];
end
if nargin < 6
    ratio = 6;
end

if isempty(sigma1)
    sigma1 = sqrt(delay1);
end
if isempty(sigma2)
    sigma2 = sqrt(delay2);
end
if isempty(t)
    t = 0: 0.05: (delay2 + 5*sigma2);   % default as in FSL
end

%%
a1 = (delay1/sigma1)^2;
b1 = (sigma1^2)/delay1;

a2 = (delay2/sigma2)^2;
b2 = (sigma2^2)/delay2;

hrf = gampdf(t, a1, b1) - gampdf(t, a2, b2)/ratio;
kernel.kernel = hrf;
kernel.t = t;
kernel.dt = t(2) - t(1);
end % function
