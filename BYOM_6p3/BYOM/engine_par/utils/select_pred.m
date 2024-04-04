function [conf_type,fname_prof,par,Tend,Cw] = select_pred(opt_in)

% Usage: [conf_type,fname_prof,par,Tend,Cw] = select_pred(opt_in)
% 
% This is a small helper function to allow the user to select MAT files for
% model curves with CIs, and exposure profiles for LPx/EPx predictions. 
% opt=1 to only load a MAT file
% opt=2 to load MAT and an exposure profile as txt file.
% opt=3 to load MAT and the data set saved as separate MAT file 
% opt=[* 1] to also load the parameter structure par from file
% opt=[* 2] to also load the parameter structure par from file, and to use
%       the parameters for the second calibration data set (with glo.names_sep)
%       NOTE: watch out as the CIs would still be made with the first set!
% opt = [* * 1] to load the complete global settings with 'glo' and 'X0mat'
% opt = [* * 2] load complete global settings, but keep pre-defined exposure scenarios
%       (if no exposure scenarios are pre-defined, this option will remove
%       any scenarios from the MAT file).
% opt = [* * * x>0] to also create an exposure profile with make_scen, using type x
% 
% Outputs are:
% conf_type  the type of sample loaded (to be used to set opt_conf.type)
% fname_prof containing the filename for the exposure profile (to be used as input for the LPx/EPx predictions)
% par        parameter structure from saved set
% Tend       gives the last time point of the exposure profile (to define glo.t)
% Cw         the exposure profile itself (to be used with make_scen for test design)
% 
% Note that this function sets glo.mat_nm, which triggers load_rnd.m to
% load the selected mat file, rather than the one based on the name of the
% main script (in glo.basenm). Also note that this function sets glo.sel,
% glo.moa and glo.feedb, if these settings are included in the file name.
% This is used by some GUTS and DEBtox analyses. Note that this function
% sets glo.Tbp that specifes the brood-pouch delay used to generate the MAT
% file (only for specific DEBtox analyses). If par is requested as output,
% this function will also set glo.names_sep if it exists.
% 
% Author     : Tjalling Jager 
% Date       : June 2022
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2022, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo DATA W X0mat

% use Matlab GUI to select MAT file for predictions
[fname_sample,filepath] = uigetfile('*.mat','Select a saved MAT file for predictions','MultiSelect','off'); 
if ~iscell(fname_sample) && numel(fname_sample) == 1 && fname_sample(1) == 0 % if cancel is pressed ...
    disp('Cancel is pressed, so no MAT file has been loaded.')
    return % simply stop
end
if ~strcmp(filepath(1:end-1),pwd)
    error('The saved MAT file must be located in the present working directory')
end

opt = zeros(1,4); % zeros will remain when opt_in is short
opt(1:length(opt_in)) = opt_in; % replace opt with input values

listOfVariables = who('-file',fname_sample); % read names of stored variables in mat file

REM    = glo;
useREM = 0;
if opt(3) == 2 && isfield(glo,'int_scen') && ~isempty(glo.int_scen) % scenarios have already been defined, and we want to keep them!
    useREM = 1;
end
 
warning('off','backtrace')
if ismember('GLO', listOfVariables)
    load(fname_sample,'GLO') % read GLO if it is saved in the mat file!
    if opt(3) > 0
        glo      = GLO; % use the complete save global settings!
    end
elseif opt(3) > 0
    error('Settings structure glo is not available in the MAT file.')
end
if ismember('X0mat', listOfVariables)
    load(fname_sample,'X0mat') % read X0mat if it is saved in the mat file!
elseif opt(3) > 0
    warning('Scenario identifiers/initial value matrix X0mat is not available in the MAT file. Make sure it is defined in the script.')
end
warning('on','backtrace')

if useREM == 1 % scenarios have already been defined, and we want to keep them!
    make_scen(-5,-1); % remove all spline info for exposure profiles (just to be on the safe side)
    disp('   exposure scenarios loaded from MAT are NOT used (pre-defined ones are kept)')
    glo.int_scen = REM.int_scen;
    glo.int_coll = REM.int_coll;
    glo.int_type = REM.int_type;
    if isfield(REM,'LabelTable') % if labels were defined, keep them as well
        glo.LabelTable = REM.LabelTable;
    end
    glo = rmfield(glo,{'wts','wts2','var'}); % remove specific fields from glo that relate to the data set
    % This keeps sameres defined as in the calibration!
end
if isfield(glo,'t')
    glo = rmfield(glo,'t'); % always remove time vector (will be re-created by prelim_checks)
end
if isfield(glo,'zvd')
    glo = rmfield(glo,'zvd'); % always remove zero-var data (will be re-created by prelim_checks)
end
% But always keep some settings as already defined
glo.basenm     = REM.basenm; 
glo.saveplt    = REM.saveplt; 

if isfield(REM,'diary') % if it was already defined ...
    glo.diary = REM.diary; % use the pre-defined one
end
if isfield(glo,'diary')
    diary (glo.diary) % collect output in the diary "results.out"
else
    diary ('results.out') % collect output in the diary "results.out"
end

disp('=================================================================================')
disp(['File selected for sample: ',fname_sample])
type = fname_sample(end-5:end-4); 
if strcmp(type,'MC')
    conf_type = 1; % Bayes slice sampler results
    disp('   Bayesian slice sample')
elseif strcmp(type,'LR')
    conf_type = 2; % likelihood region results
    disp('   likelihood-region sample')
elseif strcmp(type,'PS')
    conf_type = 3; % parspace explorer results
    disp('   parameter-space explorer sample')
elseif strcmp(type,'LP') 
    disp('   profile-likelihood results (does not include a sample!)')
    if exist([fname_sample(1:end-7),'_LR.mat'],'file') == 2
        conf_type = 2; % likelihood region results
        disp('   the similarly named likelihood-region file is available and can be used.')
    end
else
    error('The selected MAT file is not of known type.')
end

glo.mat_nm = fname_sample(1:end-7); % remove the last part (e.g., _PS.mat) of the filename so that the basename remains
% Note: setting this global implies that load_rnd will load this filename,
% rather than the one based on the script that generated it!

% We need to catch the cases where we do a DEB-TKTD analysis with an old
% MAT file, one that saved Tbp rather than the entire glo. Also when we
% don't read the entire glo, we still need to take glo.Tbp from there?
% Maybe not ... if you don't use the full glo from the MAT file you have to
% think.
if opt(3) == 0
    if ~ismember('GLO',listOfVariables) && ismember('Tbp',listOfVariables) % returns true if Tbp is stored
        load(fname_sample,'Tbp') % read glo.Tbp if it is saved in the mat file!
        glo.Tbp = Tbp;
        disp(['   Loaded from saved MAT file, brood pouch delay  glo.Tbp = ',num2str(glo.Tbp)])
%     elseif isfield(GLO,'Tbp')
%         glo.Tbp = GLO.Tbp;
%         disp(['   Loaded from saved MAT file, brood pouch delay  glo.Tbp = ',num2str(glo.Tbp)])
    end
else
    glo_entries = {'Tbp','FBV','KRV','kap','yP','Lm_ref','len','E0','delM','dV','Tref','T','sel'};
    if opt(3) == 1
        disp('   Loaded from saved MAT file, complete settings structure glo')
        disp(['   Including ',num2str(length(glo.int_scen)),' defined exposure scenarios, and, a.o.:'])
    elseif opt(3) == 2
        disp('   Loaded from saved MAT file, settings structure glo without exposure scenarios, a.o.:')
    end
    fprintf('     ')
    i_entries = 0;
    for i = 1:length(glo_entries)
        if isfield(glo,glo_entries{i})
            fprintf('glo.%s = %g ; ',glo_entries{i},glo.(glo_entries{i}))
            i_entries = i_entries + 1;
            if i_entries == 4
                fprintf('\n     ')
                i_entries = 0;
            end
        end
    end
    fprintf('\n')
    if isfield(glo,'stdDEB_start')
        disp('   Structure glo includes the starting values for stdDEB-TKTD analysis')
    end
end

% Some GUTS analyses will save the death mechanism in the file name. The
% next section looks for the text '_sel' in the mat file name, and uses it
% to define glo.sel.
k = strfind(fname_sample,'_sel'); % find where _sel is in the filename
if ~isempty(k) % if it is not empty ...
    glo.sel = str2double(fname_sample(k+4)); % this is the death mechanism belonging to the saved file
    disp(['   Loaded from saved MAT file, death mechanism    glo.sel = ',num2str(glo.sel)])
end

% Some DEBtox analyses will save the MoA and the feedback configuration in
% the file name. This looks for the text '_moa' and '_feedb' to define
% glo.moa and glo.feedb.
ku = strfind(fname_sample,'_');
k1 = strfind(fname_sample,'_moa');
if ~isempty(k1) % if it is not empty ...
    moa_tmp = fname_sample(k1+4:ku(find(ku>k1,1,'first'))-1);
    % This reads from _moa up to the next underscore!
    for i = 1:length(moa_tmp)
        glo.moa(i) = str2double(moa_tmp(i));
    end
    disp(['   Loaded from saved MAT file, mode of action     glo.moa = ',num2str(glo.moa)])
end

k2 = strfind(fname_sample,'_feedb');
if ~isempty(k2) % if it is not empty ...
    feedb_tmp = fname_sample(k2+6:ku(find(ku>k2,1,'first'))-1);
    % This reads from _feedb up to the next underscore!
    for i = 1:length(feedb_tmp)
        glo.feedb(i) = str2double(feedb_tmp(i));
    end
    disp(['   Loaded from saved MAT file, feedback config. glo.feedb = ',num2str(glo.feedb)])
end    

if opt(2) > 0 % then we also need to return par!
    if ismember('par', listOfVariables) % returns true if par is stored
        load(fname_sample,'par') % read par if it is saved in the mat file!
        disp('   Loaded from saved MAT file, parameter structure par')
    else
        error('No parameter set is present in the saved MAT file!')
    end
    
    % Some analyses contain parameters that can differ between data sets.
    % If that is the case, we need to make sure that the glo.names_sep is
    % set for the correct parameters. Otherwise, par makes no sense.
    names_sep = {};
    if ismember('names_sep',listOfVariables) % returns true if names_sep is stored
        load(fname_sample,'names_sep') % read glo.names_sep if it is saved in the mat file!
        glo.names_sep = names_sep;
    elseif exist('GLO') && isfield(GLO,'names_sep')
        names_sep = GLO.names_sep;
        glo.names_sep = names_sep;
    end
    if ~isempty(names_sep)
        disp(['   Loaded from saved MAT file, extension namesep for ',num2str(length(glo.names_sep)),' parameters'])
        if opt(2) > 1 % if the request is for a data set that is not the first ...
            disp(['   Parameters used in par are for calibration data set: ',num2str(opt(2))])
            disp('   Use with care as the CI would still be generated with the first parameter set!')
            disp('   (use option opt_ecx.id_sel instead for ECx/EPx calculations with multiple parameter sets,')
            disp('   or opt_conf.use_par_out for other cases where CIs are needed.')
            % use the parameters for data set in opt(2)
            for i_sep = 1:length(names_sep) % run through extra parameter names for separate sets
                par.(names_sep{i_sep}) = par.([names_sep{i_sep},num2str(opt(2)-1)]); % copy extra parameter level to par
            end
        end
    end
else % make sure the outputs are defined
    par = [];
end

if opt(1) == 2 % we also need to load an exposure profile
    
    % use Matlab GUI to load exposure profile from text file for predictions
    [fname_prof,filepath] = uigetfile('*.txt','Select a text file with exposure profile for predictions','MultiSelect','off');
    if ~iscell(fname_prof) && numel(fname_prof) == 1 && fname_prof(1) == 0 % if cancel is pressed ...
        disp('Cancel is pressed, so no exposure profile file has been loaded.')
        return % simply stop
    end
    
    disp(['File selected for exposure profile: ',fname_prof])
    % if the profile is not located in the working directory, return the entire path
    if ~strcmp(filepath(1:end-1),pwd)
        disp(['   in directory: ',filepath])
        fname_prof = [filepath,fname_prof];
    end
    
    if isfield(glo,'int_scen') % scenarios have already been defined
        make_scen(-5,-1); % remove all spline info for exposure profiles (just to be on the safe side)
        disp('   previously defined/loaded exposure scenarios are deleted!')
    end
    
    % Load exposure profile from file to use with linear interpolation.
    % This is needed here to define a time vector glo.t.
    Cw   = load(fname_prof);
    disp(['   exposure matrix with ',num2str(size(Cw,1)),' time points and ',num2str(size(Cw,2)-1), ' scenarios'])
    Tend = Cw(end,1); % end of the profile (for defining glo.t)
    if size(Cw,2) == 2 && Cw(2,1) ~= 0 % if the second time element is a zero, we have a file defined with identifiers
        % so here we have an exposure profile without header as e.g., a FOCUS profile
        Cw = [1 1;Cw]; % this add a first row with scenario identifier
        disp('   the exposure scenario gets the identifier 1')
    end
    
    glo.t = linspace(0,Tend,max(100,2*size(Cw,1))); % also predefine glo.t
    disp(['   time vector glo.t predefined as 0-',num2str(Tend),' with ',num2str(length(glo.t)),' points'])
    if opt(4) > 0 % then we'll also define the exposure scenario here
        make_scen(opt(4),Cw); % create exposure profile with opt(4) as type
        disp('   exposure scenario defined in global (no need to call make_scen again)')
    end
    
else % make sure the outputs are defined
    fname_prof = [];
    Cw         = [];
    Tend       = [];
end

if opt(1) == 3 % we also need to load a data set
    fname_data = glo.mat_nm; % start with the mat_nm, which has remove the _PS.mat already
    if ~isempty(k) % if it is not empty we have a GUTS analysis with _sel in the name
        fname_data = [fname_data(1:k-1),'_DATA']; % remove the last part
    elseif ~isempty(k1) % if it is not empty we have a DEBtox analysis with _moa in the name
        fname_data = [fname_data(1:k1-1),'_DATA']; % remove the last part
    else
        fname_data = [fname_data,'_DATA'];
    end
    if exist([fname_data,'.mat'],'file') ~= 2
        error('No data set with the right file name was found.')
    else
        listOfVariables = who('-file',fname_data); % read names of stored variables in mat file
        load(fname_data,'DATA','W') % read DATA and W if saved in the mat file!
        disp('   Loaded from saved MAT file, data set')
        if ismember('DATAx', listOfVariables)
            global DATAx Wx
            load(fname_data,'DATAx','Wx') % read DATAx and Wx if saved in the mat file!
            disp('   Loaded from saved MAT file, extra data set')
        end
    end
end

disp('=================================================================================')
diary off  % close the diary that collected the output to screen
