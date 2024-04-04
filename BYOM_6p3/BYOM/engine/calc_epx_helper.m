function [crit,Xout] = calc_epx_helper(MF,calc_int,t,par,X0v,glo,Xctrl,locX,Feff,WRAP2)

% Usage:  [crit,Xout] = calc_epx_helper(MF,calc_int,t,par,X0v,glo,Xctrl,locX,Feff,WRAP2)
% 
% Small helper function to calculate results needed for effect window and
% EPx calculations. This works to output effects at the end of <t>,
% relative to the control, but also a criterion to use with <fzero> to find
% exact EPx values. This is code that used to be in sub-functions of the
% effect window and EPx functions.
% 
% Inputs:
% <MF>       the multiplication factor to use on the exposure profile
% <calc_int> integrate survival and repro into 1) RGR, or 2) survival-corrected repro (experimental!)
% <t>        the time vector for calculating the response of the traits
% <par>      the optimised parameter set
% <X0v>      vector with scenario and initial values for the state variables
% <glo>      structure with information that is generally global
% <Xctrl>    value for the control trait (as EPx is relative to control)
% <locX>     location of trait(s) of interest in Xout
% <Feff>     the fraction effect (x/100 in EPx)
% <WRAP2>    a wrapper containing various variables needed for the intrinsic rate
% 
% Author     : Tjalling Jager 
% Date       : August 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2022, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

crit   = NaN; % make sure this output is defined
glo.MF = MF; % modify global with multiplication factor
X = call_deri(t,par,X0v,glo); % use call_deri.m to provide the model output

%% If we need to return single traits ...

if calc_int == 0 % need to return one or more traits, relative to control
    
    Xout = X(end,locX)./Xctrl; % return final values of selected trait(s), relative to control response
    
    if ~isempty(Feff)
        if length(Xout)>1 % this should not be needed since this function is only called by other functions
            error('Use fzero only for 1 endpoint!')
        end
        crit = Xout-(1-Feff); % zero when end value is x*100% effect
    end
    
    return % return to calling function
end

%% If we need to return integrated traits ...

Rc = X(:,glo.locR); % take cumulative repro from the output
S  = X(:,glo.locS); % take survival from the output

if isempty(S)
    S = ones(size(Rc)); % just assume no deaths
end

if calc_int == 1 % intrinsic rate of increase
    
    % unwrap the wrapped variables
    t  = WRAP2.t;
    t2 = WRAP2.t2;
    Th = WRAP2.Th;
    rgr_init = WRAP2.rgr_init;
    
    rgr = NaN;
    Rp  = diff(Rc)./(diff(t)); % estimated repro rate, as difference (on time vector t2)
    Sp  = mean([S(1:end-1) S(2:end)],2); % mean survival probability in an interval delta t (on time vector t2)
    if any(Rp>0) % only if there is repro (otherwise, the initial NaN remains)
        % create an anonymous function with the current values for the traits
        crit = @(rgr) trapz(t2,Rp .* Sp .* exp(-rgr*(t2+Th))) -1; % trapezium rule integration
        % extract 1 as we need to find where the integrated thing is 1
        if crit(0) <= 0
            rgr = 0; % don't calculate negative rgr's
        else
            rgr = fzero(crit,rgr_init); % use fzero to find the RGR
        end
    end
    Xout = rgr/Xctrl;
    
else % calculate survival-corrected reproduction as endpoint.
    
    Rp    = diff(Rc); % estimated repro in each interval delta t
    Sp    = mean([S(1:end-1) S(2:end)],2); % mean survival probability in each interval delta t
    Rcorr = sum(Rp .* Sp); % summed survival-corrected repro in each interval
    % we only need the last value for the cumulated reproduction, so sum is enough
    Xout  = Rcorr/Xctrl;
    
end

if ~isempty(Feff)
    crit = Xout-(1-Feff); % zero when end value is x*100% effect
end
