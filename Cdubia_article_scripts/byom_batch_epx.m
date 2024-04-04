%% %% BYOM, ERA-special with simple-compound model: byom_batch_epx.m
%
% *Table of contents*

%% About
% * Author     : Tjalling Jager
% * Date       : December 2022
% * Web support: <http://www.debtox.info/byom.html>
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
% (see <http://www.debtox.info/book_debkiss.html>) provides a description
% of the model, as well as the publication of Jager in Ecological
% Modelling: <https://doi.org/10.1016/j.ecolmodel.2019.108904>.
%
% *This script:* This script demonstrates the use of the parameter-space
% explorer from the openGUTS project (see <http://www.openguts.info/>) in
% making batch-window predictions for FOCUS profiles. It runs through all
% selected profiles to calculate EPx for each window. Only for the lowest
% EPx, a CI will be calculated (unless CIs are suppressed). This script
% works generically, as all species- and compound-specific information is
% loaded from the MAT file. Furthermore, it is now even largely model
% independent, as glo and X0mat are also saved during the calibration, and
% loaded here. All DEB-based TKTD analyses should work with this script.
% 
%  Copyright (c) 2012-2023, Tjalling Jager, all rights reserved.
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
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:

% No data for simulations!

%% Call the Matlab GUI open-file element to load profile and MAT file

[conf_type,~,par] = select_pred([1 1 1]); 
par.a = [1 0 1e-3 1e3 1];
% Note that this function, with these settings, loads a pre-saved MAT file.
% It loads the parameters, as well as all settings used to generate the MAT
% file (in glo) apart from the exposure profile (which we'll redefine
% anyway). This only works when the MAT file was generated by BYOM v.6.0
% BETA 7 or newer as older versions won't save glo. When using an old MAT
% file (pre BYOM v6), make sure to set glo.Tbp to the correct value if it
% needs to be >0.

% The parameter structure par is loaded from the selected MAT file. Even
% though it does not need to be used to call other functions, it must be
% specified for the code to run error free.

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios. All identifiers and initial values, as part of
% the matrix X0mat, are also loaded by select_pred already. Feel free to
% overwrite them with other values if needed.
% 
% NOTE: the glo.locD etc. have all been loaded from the MAT file by select_pred..

%% Initial values for the model parameters
% All parameters have been loaded from file using select_pred. All global
% parameters, as part of the structure glo, are also loaded by select_pred
% already. Feel free to overwrite them with other values.
  
% glo.len    = 2;       % switch to fit physical length (0=wwt, 1=phys. length, 2=phys. length, no shrinking) (used in call_deri.m)

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

glo.t = linspace(0,10,100); % need to define time vector as we have no data
% Note: the setting of this global has no effect on the results of the
% analyses below, since the time vector is redefined for calculating ECx
% and EPx anyway.

% All plot labels have been loaded from file using select_pred.

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
% well. Also note that glo.stiff and glo.break_time are loaded from the MAT
% file as well, but overwritten below. For FOCUS exposure scenarios, we
% best use glo.break_time=0, and possibly a different solver than for
% constant exposure.
 
% -------------------------------------------------------------------------
% Configurations for the ODE solver
glo.stiff = [0 3]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Second argument is for default sloppy (0), normally tight (1), tighter
% (2), or very tight (3) tolerances. Use 1 for quick analyses, but check
% with 3 to see if there is a difference! Especially for time-varying
% exposure, there can be large differences between the settings!
glo.break_time = 0; % break time vector up for ODE solver (1) or don't (0)
% NOTE: for FOCUS profiles, breaking the time vector does not help accuracy
% much, but it does slow down the calculations a lot.
% -------------------------------------------------------------------------

% No need for optimisation or standard plotting here.

%% Analyse exposure profiles to check safety margins
% The functions below analyse the impacts of a series of exposure profiles
% (e.g., from FOCUS). First, a batch-mode 'moving time-window' analysis is
% demonstrated, with a fixed window, calculating explicit EPx for each
% window. Only for the lowest EPx in each profile, a CI is calculated (when
% opt_conf.type > 0).
% 
% Note that the functions calc_epx and calc_epx_window use the options
% structure opt_ecx. This option structure, by default, makes sure that the
% parameter called _hb_ is set to zero in the analysis (both in the
% best-fit parameter set and in the sample).

opt_conf.type     = conf_type; % make intervals from 1) slice sampler, 2) likelihood region, 3) parspace explorer
opt_conf.lim_set  = 2; % use limited set of n_lim points (1) or outer hull (2, not for Bayes) to create CIs
opt_conf.n_lim    = 100; % size of limited set (likelihood-region and parspace only)
opt_ecx.statsup   = [glo.locL glo.locS]; % states to suppress from the calculations (e.g., glo.locL)
opt_tktd.statsup   = [glo.locL glo.locS]; % states to suppress from the calculations (e.g., glo.locL)

opt_ecx.id_sel    = [0 1 0]; % scenario to use from X0mat, scenario identifier, flag for ECx/EPx to use scenarios rather than concentrations
% opt_ecx.id_sel    = [0 300 1]; % scenario to use from X0mat, scenario identifier, flag for ECx/EPx to use scenarios rather than concentrations

opt_ecx.rob_win   = 0; % set to 1 to use robust EPx calculation rather than with fzero
% Note: robust EPx calculation calculates EPx at the given steps only and
% then interpolates. This is not very precise, but with carefully selected
% steps it will be the lowest. In rather extreme cases, there will be more
% than one EPx for a specific window and a specific endpoint (when there
% are effects on growth, and feedbacks affecting k_d; a warning will be
% given). The 'regular' method may end up in either EPx (or produces an
% error), while 'robust' has a far better chance to yield the lowest.
% Robust is only used for EPx calculation, not for effect windows with
% fixed MFs.

% % Force EPx calculations to use f=1
% par.f(1) = 1;
% opt_conf.use_par_out = 1; % set to 1 to use par as entered into the function for CIs, rather than from saved set

opt_ecx.batch_plt = 1; % when set to 1 create plots when using calc_epx_window_batch
opt_ecx.showprof  = 1; % set to 1 to show exposure profile in plots at top row (calc_epx_window)
opt_ecx.start_neg = 1; % set to 1 to start the moving window at minus window width
opt_ecx.prune_win = 1; % set to 1 to prune the windows to keep the interesting ones
% Note: this setting screens all windows to find the window with highest
% minimum concentration. Any window whose maximum is lower than this value
% can be ignored. This is absolutely not guaranteed to work when there
% are effects on growth, and feedbacks affecting k_d (a warning will be
% given when trying pruning).

Twin = 9; % length of the time window (one element)
opt_ecx.Tstep    = 1; % stepsize or resolution of the time window (default 1 day)
opt_ecx.calc_int = 0; % integrate survival and repro into 1) RGR, or 2) survival-corrected repro (experimental!)

% opt_ecx.saveall = 1;   % set to 1 to save all output (EPx for each element of the sample) in separate folderr
% batch calculation of EPx
opt_ecx.Feff = [0.10]; % effect levels (>0 en <1), x/100 in ECx/EPx

%calc_epx_window_batch(par,Twin,opt_ecx,opt_conf);


% in a loop analyse all the profiles
files = dir(fullfile("profiles",'/*.txt'));  % get recursively in all the directories the profile files
for i = 1:length(files)
    [MinColl,MinCI,ind_traits]=calc_epx_window(par,fullfile(files(i).folder,files(i).name),Twin,opt_ecx,opt_conf);
    [~,ind_min] = min(MinColl{opt_ecx.Feff==0.1}(:,2)); % find where lowest EP10 is 
    Tstart = MinColl{opt_ecx.Feff==0.1}(ind_min,3); % start time for window with lowest EP10
    Trng = [Tstart, Tstart + Twin];
    calc_epx(par,fullfile(files(i).folder,files(i).name), ...
            Trng,opt_ecx,opt_conf,opt_tktd)
end