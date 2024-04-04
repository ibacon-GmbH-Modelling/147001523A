function save_rnd_txt

% This function loads a saved MAT file (currently only when saved by (a
% recent version of) the parameter-space explorer. It then saves the outer
% rim to a text file, with all parameters on normal scale. It only saves
% the fitted parameters. This file can then be used to generate CIs in
% other software implementations. This function works outside of BYOM as
% well; it does not need anything from the engine folder or globals etc.
%
% Author     : Tjalling Jager 
% Date       : November 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2022, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

% use Matlab GUI to select MAT file for predictions
[fname_sample,filepath] = uigetfile('*_PS.mat','Select a saved parspace MAT file for predictions','MultiSelect','off'); 
if ~iscell(fname_sample) && numel(fname_sample) == 1 && fname_sample(1) == 0 % if cancel is pressed ...
    return % simply stop
end

% This now only works for mat files generated with parspace explorer!
load([filepath,fname_sample],'par','pmat','coll_all')

crit_add = 0.3764; % small addition to chi2 criterion for inner rim (to make sure coverage is adequate)
                   % 0.3764 is in line with the value used in the openGUTS software

% for the sample of the likelihood-based joint confidence region, we make a
% small addition to chi2 criterion for the inner rim (to make sure coverage
% is adequate, since we use a discrete, limited, sample)
chicrit   = 3.8415+crit_add; % critical value for 95% confidence at df=1, with an addition
chicrit2  = 3.8415-2*crit_add; % with a subtraction ...

% this loads the sets from the conf region in coll_all, and the parameter
% matrix pmat. Note that coll_all is already sorted based on likelihood.
%
% last column in coll_all is log-likelihood itself

% Modern BYOM version that also saves par
% In future versions, this will likely remain as the only
% method to load a sample.

mll = coll_all(1,end); % extract maximum likelihood value from saved set (<coll_all> is sorted, so it is in the first row)
rnd = coll_all; % rename coll_all to rnd
clear coll_all; % clear coll_all to save memory

rnd(:,end) = 2 * (rnd(:,end)-mll);       % translate into 2 times loglik ratio
rnd        = rnd(rnd(:,end)<chicrit,:);  % keep the sets in inner rim (df=1)
rnd        = rnd(rnd(:,end)>chicrit2,:); % keep the sets in outer part of inner rim (df=1)
rnd        = rnd(:,1:end-1);             % and remove last column

% In <rnd> the parameters that are fitted on log-scale are still on
% log-scale. We'll put them on normal scale below.
pmat       = pmat(pmat(:,2)==1,:); % remove non-fitted parameters
ind_logfit = pmat(:,5) == 0; % these are the parameters that are fitted on log scale

if sum(ind_logfit) > 0
    rnd(:,ind_logfit) = 10.^(rnd(:,ind_logfit)); % put on normal scale
end

% save as text file
save([filepath,fname_sample(1:end-7),'_outer_rim.txt'],'rnd','-ascii','-tabs') % save text file in same directory as the _PS.mat!

% Also save a text file with the best values and the fixed parameters
names  = fieldnames(par); % extract field names from par
fileID = fopen([filepath,fname_sample(1:end-7),'_best_vals.txt'],'w');
for i = 1:length(names)
    fprintf(fileID,'%s %0.5g %1.0f\n',names{i},par.(names{i})(1:2));
end
fclose(fileID);    