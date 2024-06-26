%% BYOM, ERA-special-Daphnia with simple-compound model: byom_debtox_ceriodaphnia.m
%
% *Table of contents*

%% About
% * Author     : Tjalling Jager
% * Date       : November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_debtox2019.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* simple DEBtox model for toxicants, based on DEBkiss and
% formulated in compound parameters. The model includes flexible modules
% for toxicokinetics/damage dynamics and toxic effects. The DEBkiss e-book
% (see <http://www.debtox.info/book_debkiss.html>) provides a partial
% description of the model; the publication of Jager in Ecological
% Modelling contains the full details:
% <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
%
% *This script:* The water flea _Daphnia magna_ exposed to fluoranthene;
% data set from Jager et al (2010),
% <http://dx.doi.org/10.1007/s10646-009-0417-z>. Published as case study in
% Jager & Zimmer (2012). Simultaneous fit on growth, reproduction and
% survival. Parameter estimates differ slightly from the ones in Jager &
% Zimmer due to several model changes and choices in the analysis.
% Furthermore, here the raw data on individual females is used, such that
% the special functions for censoring the reproduction data can be used.
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(1) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 1; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% NOTE: time MUST be entered in DAYS for the estimation of starting values
% to provide proper search ranges! Controls must use identifiers 0 for the
% true control and 0.1 for the solvent control. If you want to use another
% identifier for the solvent, change parameter id_solvent in
% automatic_runs.

TR   = 0.5; % transformations for continuous data
opt  = 1; % select an option (opt=1 is recommended for Daphnia)

% Options to deal with repro data set:
% 0) Check whether we can use a single intermoult period for the entire
%    data set. Screen output will show mean intermoult times and brood 
%    sizes across the replicates, as function of brood number and treatment.
% 1) Cumulate reproduction, but remove the time points with zero
%    reproduction. Good for clutch-wise reproduction.
% 2) Cumulate reproduction, but don't remove zeros. Good for continuous
%    reproduction or when animals are not followed individually.
% 3) Shift neonate release back to the previous moult. When this option is
%    used, don't shift the model predictions with <glo.Tbp>: the data now
%    represent egg production rather than neonate release.

% Put the position of the various states in globals, to make sure that the
% correct one is selected for extra things (e.g., for plotting in
% plot_tktd, in call_deri for accommodating 'no shrinking', for population
% growth rate).
glo.locD = 1; % location of scaled damage in the state variable list
glo.locL = 2; % location of body size in the state variable list
glo.locR = 3; % location of cumulative reproduction in the state variable list
glo.locS = 4; % location of survival probability in the state variable list

fnames = {'data_AZT_Cdubia_calibration'};
read_datafiles(fnames,TR,opt); % helper function defines DATA, W and glo.LabelTable

% Note: optionally, add some info to the MAT filename. The MAT filename
% will already include MoA and feedbacks, but if you want to try other
% things as well (changing opt, calibrating on the validation data, etc)
% it can be helpful to change the name to use this script but not 
% overwrite previous MAT files.

% glo.basenm  = [mfilename,'_CAL1']; % remember the filename for THIS file for the plots
save([glo.basenm,'_DATA'],'DATA','W') % save MAT file with data set 
% Saving the data set is handy to allow for a simple reconstruction of the
% calibrations, without needing to define the data again, in the same way.
% Note that if you use DATAx and Wx, you need to save them as well!

switch opt
    case 0 % check intermoult duration (and mean brood size)
        return % we need to stop to check the results (no data are created)
    case 3 % then we have shifted the data set, so no shifting of model needed
        glo.Tbp = 0; % time that eggs spend in brood pouch (model output for repro will be shifted)
    otherwise % we need to shift the model output
        %glo.Tbp = 1; % time that eggs spend in brood pouch (model output for repro will be shifted)
        glo.Tbp = 0; % in some cases, effect may be on the eggs in the pouch, so no shift proposed
end

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat      = [glo.LabelTable.Scenario]'; % the scenarios (here identifiers) 
X0mat(2,:) = 0; % initial values state 1 (scaled damage)
X0mat(3,:) = 0; % initial values state 2 (body length, initial value overwritten by L0)
X0mat(4,:) = 0; % initial values state 3 (cumulative reproduction)
X0mat(5,:) = 1; % initial values state 4 (survival probability)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 
  
% global parameters as part of the structure glo
glo.FBV    = 0.02;    % dry weight egg as fraction of structural body weight (-) (for losses with repro; approx. for Daphnia magna)
glo.KRV    = 1;       % part. coeff. repro buffer and structure (kg/kg) (for losses with reproduction)
glo.kap    = 0.8;     % approximation for kappa (for starvation response)
glo.yP     = 0.8*0.8; % product of yVA and yAV (for starvation response)
glo.Lm_ref = 1;       % reference max length for scaling rate constants
glo.len    = 2;       % switch to fit length 1) with shrinking, 2) without shrinking (used in call_deri.m)
% NOTE: the settings above are species specific! Make sure to use the same
% settings for validation and prediction as for calibration! 
% 
% NOTE: For arthropods, one would generally want to fit the model without
% shrinking (since the animals won't shrink in length).

% mean length at start taken as mean of the true controls at t=0, and fixed below
mL0 = mean(DATA{1,2}(2,1+find(DATA{1,2}(1,2:end)==0)));

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.L0   = [mL0   0 0.2  0.4 1]; % body length at start experiment (mm)
par.Lp   = [0.6646  1 0.4    1  1]; % body length at puberty (mm)
par.Lm   = [0.8842     1 0.5    2  1]; % maximum body length (mm)
par.rB   = [0.2908   1 0.1    1  1]; % von Bertalanffy growth rate constant (1/d)
par.Rm   = [9.482    1   3   20  1]; % maximum reproduction rate (#/d)
par.f    = [1     0   0    2  1]; % scaled functional response
par.hb   = [0.001  0 1e-3 0.07 1]; % background hazard rate (d-1)
% - Note 1: it does not matter whether length measures are entered as actual
% length or as volumetric length (as long as the same measure is used
% consistently). 
% - Note 2: hb is fitted on log-scale. This is especially helpful for
% fit_tox(1)=-2 when hb is fitted along with the other parameters. The
% simplex fitting has trouble when one parameter is much smaller than
% others. 
% - Note 3: the min-max ranges in par are appropriate for Daphnia magna. For
% other species, these ranges and starting values need to be modified.

glo.names_sep = {}; % no parameters can differ between data sets
% glo.names_sep = {'f';'L0'}; % names of parameters that can differ between data sets
% par.L01   = [0.9 1 0.5  1.5  1]; % body length at start experiment (mm)
% par.f1    = [0.9 1 0    2    1]; % scaled functional response

% Note: using separate parameters for separate data sets requires using
% specific identifiers, and using exposure scenarios with make_scen (also
% for constant exposure and for the controls)!

% extra parameters for special situations
par.Lf   = [0 0 0 1e6 1]; % body length at half-saturation feeding (mm)
par.Lj   = [0 0 0 1e6 1];  % body length at end acceleration (mm)
par.Tlag = [0 0 0 1e6 1];  % lag time for start development

ind_tox = length(fieldnames(par))+1; % index where tox parameters start
% the parameters below this line are all treated as toxicity parameters!

% When using the parameter-space explorer, start values and ranges are not
% needed (filled later by startgrid_debtox); only the fit/fix mark is
% relevant here. But make sure that the value in the first column is within
% the bounds (and not zero for log-scale parameters).
par.kd   = [0.08   1 0.01  10 0]; % dominant rate constant (d-1)
par.zb   = [0.1    1 0    1e6 1]; % effect threshold energy budget ([C])
par.bb   = [75     1 1e-6 1e6 0]; % effect strength energy-budget effect (1/[C])
par.zs   = [0   0 0    1e6 1]; % effect threshold survival ([C])
par.bs   = [1e-6    0 1e-6 1e6 0]; % effect strength survival (1/([C] d))
 
% After optimisation, copy-paste relevant lines (fitted parameters) from screen below!

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = ['scaled damage (',char(181),'g/L)'];
glo.ylab{2} = 'body length (mm)';
if isfield(glo,'Tbp') && glo.Tbp > 0
    glo.ylab{3} = ['cumul. repro. (shift ',num2str(glo.Tbp),'d)'];
else
    glo.ylab{3} = 'cumul. repro. (no shift)';
end
glo.ylab{4} = 'survival fraction (-)';

% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'conc. '; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'M']; % legend label after the 'scenario' number
% Note: these legend labels will not be used when we make a glo.LabelTable

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 
% 
% NOTE: for this package, the options useode and eventson in glo will not
% be functional: the ODE solver is always used, and the events function as
% well.
 
% -------------------------------------------------------------------------
% Configurations for the ODE solver
glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for default sloppy (0), normally tight (1), tighter
% (2), or very tight (3) tolerances. Use 1 for quick analyses, but check
% with 3 to see if there is a difference! Especially for time-varying
% exposure, there can be large differences between the settings!
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% Note: breaking the time vector is a good idea when the exposure scenario
% contains discontinuities. Don't use for continuous splines (type 1) as it
% will be much slower. For FOCUS scenarios (high time resolution), breaking
% up is not efficient and does not appear to be necessary.
% -------------------------------------------------------------------------

opt_optim.fit    = 1; % fit the parameters (1), or don't (0)
opt_plot.bw      = 0; % if set to 1, plots in black and white with different plot symbols
opt_plot.annot   = 2; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
opt_plot.repls   = 1; % set to 1 to plot replicates, 0 to plot mean responses
basenm_rem       = glo.basenm; % remember basename as automatic_runs may modify it!

% Select what to fit with fit_tox (this is a 3-element vector). NOTE: use
% identifier c=0 for regular control, and c=0.1 for solvent control. When
% entering more than one data set, use 100 and 100.1 for the controls of
% the second data set (and 101, 102 ... for the treatments), 200 and 200.1
% for the controls of the third data set. etc. 
% 
% First element of fit_tox is which part of the data set to use:
%   fit_tox(1) = -2  comparison between control and solvent control (c=0 and c=0.1)
%   fit_tox(1) = -1  control survival (c=0) only
%   fit_tox(1) = 0   controls for growth/repro (c=0) only, but not for survival
%   fit_tox(1) = 1   all treatments, but, when fitting, keep all control parameters fixed; 
%               run through all elements in MOA and FEEDB sequentially and 
%               provide a table at the end (plots are made incl. control)
% 
% Second element of fit_tox is whether to fit or only to plot:
%   fit_tox(2) = 0   don't fit; for standard optimisations, plot results for
%               parameter values in [par], for parspace optimisations, use saved mat file.
%   fit_tox(2) = 1   fit parameters
% 
% Third element is what to use as control (fitted for fit_tox(1) = -1 or 0)
% (if this element is not present, only regular control will be used)
%   fit_tox(3) = 1   use regular control only (identifier 0)
%   fit_tox(3) = 2   use solvent control only (identifier 0.1)
%   fit_tox(3) = 3   use both regular and control (identifier 0 and 0.1)
% 
% The strategy in this script is to fit hb to the control data first. Next,
% fit the basic parameters to the control data. Finally, fit the toxicity
% parameter to the complete data set (keeping basic parameters fixed. The
% code below automatically keeps the parameters fixed that need to be
% fixed. 
% 
% MOA: Mode of action of toxicant as set of switches
% [assimilation/feeding, maintenance costs (somatic and maturity), growth costs, repro costs] 
% [1 0 0 0 0]   assimilation/feeding
% [0 1 0 0 0]   costs for maintenance 
% [0 0 1 1 0]   costs for growth and reproduction
% [0 0 0 1 0]   costs for reproduction
% [0 0 0 0 1]   hazard for reproduction
% 
% FEEDB: Feedbacks to use on damage dynamics as set of switches
% [surface:volume on uptake, surface:volume on elimination, growth dilution, losses with reproduction] 
% [1 1 1 1]     all feedbacks 
% [1 1 1 0]     classic DEBtox (no losses with repro)
% [0 0 1 0]     damage that is diluted by growth
% [0 0 0 0]     damage that is not diluted by growth

% These are the MoA's and feedback configurations that will be run
% automatically when fit_tox(1) = 1. This needs to be defined here, even
% for control fits when it is not used.

%MOA   = [1 0 0 0 0;0 1 0 0 0;0 0 1 1 0;0 0 0 1 0];
FEEDB = [0 0 0 0; 1 0 0 0; 1 1 0 0; 1 1 1 0; 0 0 1 0];
MOA= [0 0 0 0 1];
% MOA   = [0 0 0 1 0; 0 0 0 0 1]; % these are costs for repro and repro hazards
% FEEDB = [0 0 0 0; 1 0 0 0;1 1 1 0]; % you may want to try more feedbacks ...
% FEEDB = allcomb([0 1],[0 1],[0 1],[0 1]); % or test ALL feedback configurations

% For fit_tox(1)~=1, the setting of MOA and FEEDB has no impact; they must
% be defined to prevent errors. For fit_tox(1)=1, setting is relevant. Note
% that if you run multiple configurations, the last one will remain in the
% memory, unless glo.moa and glo.feedb are redefined. 

% Note: when using the parspace explorer, the automatic_runs uses the start
% ranges as produced by startgrid_debtox. These are preliminary! It may be
% needed to tweak these, but that should then be done in startgrid_debtox
% or in automatic_runs. So watch out when the analysis runs into min-max
% bounds.

% ===== FITTING CONTROLS ==================================================
% Simplex fitting works fine for control data
opt_optim.type = 4; % optimisation method: 1) default simplex, 4) parspace explorer

% opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
% opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
% opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)

% % Compare controls in data set
% fit_tox = [-2 1 3];
% automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot);
% % script to run the calculations and plot, automatically
% % NOTE: I think the likelihood-ratio test is often too strict, and that
% % using both controls should be the default situation.

% % Fit hb in data set
% fit_tox = [-1 1 2]; % use solvent control only
% %par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot); 
% par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot,opt_prof); 
% % script to run the calculations and plot, automatically
% par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% Fit other control parameters (not hb) in data set
fit_tox = [0 1 2]; % use solvent control only
%par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot);
par_out = automatic_runs_debtox2019(fit_tox,par,ind_tox,[],MOA,FEEDB,opt_optim,opt_plot,opt_prof);
% script to run the calculations and plot, automatically
par = copy_par(par,par_out,1); % copy fitted parameters into par, and keep fit mark in par

% ===== FITTING TOX DATA ==================================================
% Use the parameter-space explorer for fitting the treatments. Note that,
% with fit_tox=1, the tox parameters in par are replaced (when fitted) with
% estimates based on the data set using startgrid_debtox (called in
% automatic_runs_debtox2019).
% 
% Note: you can use the rough settings to find better ranges and restart
% with refined settings. With opt_optim.ps_profs = 0, the algorithm will
% provide new search ranges on screen that can be directly copied-pasted
% into this script. You may comment out the fitting of the controls above.
% Make sure that skip_sg is then set to 1 to avoid startgrid_debtox to
% overwrite par again.

opt_optim.type     = 4; % optimisation method 1) simplex, 4 parameter-space explorer
% Note: to use saved set for parameter-space explorer, use fit_tox(2)=0!
opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots of parameter space to monitor progress
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)
% -----------------------------------------------------------------------
% THIS BLOCK: SETTINGS FOR ROUGHLY GOING THROUGH MANY CONFIGURATIONS
% opt_optim.ps_rough = 2; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS (2 for extra rough)
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
% glo.stiff          = [2 1]; % for quick exploration of many options, could try stiff solver with sloppy tolerances if default setting is slow/stuck
glo.stiff          = [0 3]; % use ode45 with very strict tolerances
skip_sg            = 0; % set to 1 to skip startgrid completely (use ranges in <par> structure)
% -----------------------------------------------------------------------
% % THIS BLOCK: SETTINGS TO REFINE FROM RANGES COPIES FROM SCREEN INTO THIS SCRIPT
% opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
% glo.stiff          = [0 3]; % use ode45 with very strict tolerances
% skip_sg            = 1; % set to 1 to skip startgrid completely (use ranges in <par> structure)
% -----------------------------------------------------------------------

fit_tox = [1 1 3]; % use both controls, and fit tox data
[par_out,best_MoaFb] = automatic_runs_debtox2019(fit_tox,par,ind_tox,skip_sg,MOA,FEEDB,opt_optim,opt_plot);
% script to run the calculations and plot, automatically
% Note: automatic_runs will return the BEST parameter set in par_out.
disp_settings_debtox2019 % display a bit more info on the settings on screen
print_par(par_out) % print the complete best parameter vector that can be copied-pasted
% This includes the control parameters (the parspace explorer itself will
% only plot the results for the fitted parameters).

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% effects data, which are more readable when plotting with various
% intervals.

opt_conf.type    = 3; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_tktd.repls   = 1; % plot individual replicates (1) or means (0)
opt_tktd.transf  = 1; % set to 1 to calculate means and SEs including transformations
opt_tktd.obspred = 1; % plot predicted-observed plots (1) or not (0)
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
opt_tktd.sppe    = 1; % set to 1 to calculate SPPEs (relative error at end of test)
opt_tktd.statsup = 4;

% change X0mat to avoid plotting controls not used for fitting (assumes
% that regular control has identifier 0 and solvent control 0.1)
switch fit_tox(3)
    case 1 % remove solvent control
        X0mat(:,X0mat(1,:)==0.1) = [];
    case 2 % remove regular control
        X0mat(:,X0mat(1,:)==0) = [];
    case 3
        % leave all in
end

% Below some tricks to allow plotting results for selected configurations.
glo.moa    = MOA(best_MoaFb(1),:); % change global for MoA
glo.feedb  = FEEDB(best_MoaFb(2),:); % change global for feedback configuration
% glo.moa    = MOA(1,:);   % change global for MoA to first one
% glo.feedb  = FEEDB(1,:); % change global for feedback configuration to first one

glo.mat_nm = [basenm_rem,'_moa',sprintf('%d',glo.moa),'_feedb',sprintf('%d',glo.feedb)];
% The glo.mat_nm specifies the filename to read the parameters and sample
% from, as needed for plot_tktd. 
glo.basenm = basenm_rem; % return the basenm to this filename 
% The plots saved will get their filename based on the name of THIS
% script (so without the specification of MoA and feedbacks). 

plot_tktd([],opt_tktd,opt_conf); 
% leave the options opt_conf empty to suppress all CIs for these plots
% leave par_out (first input) empty to read it from saved mat file. If we
% ran more configurations, we should watch out when using par_out
% (automatic_run returns the best-fitting one).
