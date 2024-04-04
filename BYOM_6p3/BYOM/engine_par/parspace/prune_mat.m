function coll_all = prune_mat(coll_all,mll_crit,sel)

% Usage: coll_all = prune_mat(coll_all,mll_crit,sel)
% 
% This tiny piece of code goes through coll_all and removes all values that
% are really bad fits. Next, it goes through all values to see which lines
% are duplicates. The duplicates are subsequently removed.
% 
% Inputs
% <coll_all> the full matrix to prune rows from
% <mll_crit> remove entries with MLL (last column) larger than this value
% <sel>      vector with two elements, 0 or 1, for selecting culling based
%            on MLL and/or removing duplicate rows
% 
% Output
% <coll_all> the full matrix after pruning
% 
% Author     : Tjalling Jager
% Date       : August 2022
% Web support: <http://www.debtox.info/byom.html> and <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2022, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

% Prune <coll_all> to remove all values that are outside highest allowed chi2 criterion.
if sel(1) == 1 && ~isempty(mll_crit)
    coll_all = coll_all(coll_all(:,end) < mll_crit,:);
end

% Next, remove duplicates. This generally takes up to a few seconds.
% However, for really large coll_all (especially when rough=0), it can go
% up to several minutes, so don't call it when not really needed.
if sel(2) == 1
    t1 = toc;
    fprintf('Removing duplicate sets in sample ... ')
    SZ(1) = size(coll_all,1);
    i = 1;
    while i <= size(coll_all,1) % while not at the end of the matrix ...

        ind_i = ismember(coll_all,coll_all(i,:),'rows'); % check which rows are the same as this one

        if sum(ind_i)>1 % if that's more than one (there's always the one at i)
            ind_r    = find(ind_i==1); % look where they are (logicals to indices)
            ind_r(1) = [];             % keep the first one
            coll_all(ind_r,:) = [];    % but remove the duplicates
            % Note: you can also do this with logical indexing (and remove
            % the first '1'), but that does not improve speed noticably.
            % The slow thing here is ismember. It might be possible to find
            % an algorithm to do this is parallel ...

        end

        i = i + 1; % then go to the next one

    end
    SZ(2) = size(coll_all,1);
    t2 = toc;
    fprintf('%d sets removed (in %.1f sec.)\n',-diff(SZ),t2-t1)
    
end
