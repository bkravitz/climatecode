function workdemand=heatstress(WBGT,acclimatized,worktype,adder,opts)

% function workdemand=heatstress(WBGT,acclimatized,worktype,[adder],['Units',units])
% Calculates the work/rest demands at Hanford: the percentage of a work
% shift that can be worked safely (the rest being rest) for a given Wet
% Bulb Globe Temperature (WBGT), acclimatization state, and work rate.
%
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 15 November 2017
% Cleaned up 9 October 2026 (with Claude):
%   (1) vectorized: WBGT may be any array; the output has the same size
%       (previously only scalars worked);
%   (2) NaN WBGT now gives NaN (previously 0, i.e., "cannot work");
%   (3) optional name-value argument 'Units' ('F', default, or 'C'), so
%       the degrees C output of calcwbgt.m can be passed directly.  With
%       'Units','C' the adder is also in degrees C (see adder below);
%   (4) input validation (acclimatized must be 0 or 1, worktype 1-4);
%   (5) expanded help (table in degrees F and C, clothing adders,
%       references).
%   The table values and the lookup are unchanged; for scalar inputs the
%   result is identical to the 2017 version except that NaN gives NaN.
%
% Based on Hanford Climate Vulnerability Assessment and Resilience Action Plan
% MSA-1604467
% Davis et al., 2016
%
% The table below is that of Davis et al. (2016).  The acclimatized values
% equal the ACGIH Threshold Limit Values (TLV) for heat stress and the
% non-acclimatized values the ACGIH Action Limits, converted from degrees C
% to degrees F.  They are interpreted as upper limits for work schedules:
% the result is the largest work fraction whose limit is >= WBGT (+adder).
% For example, for non-acclimatized moderate work the 25% limit is 84.2 F,
% so for WBGT above 84.2 F work cannot be done safely (result 0).
%
% WBGT upper limits, degrees F (degrees C)
%
% Acclimatized (ACGIH TLV)
%   Work/rest    Light        Moderate     Heavy        Very Heavy
%   100% work    85.1 (29.5)  81.5 (27.5)  78.8 (26.0)  --
%    75% work    86.9 (30.5)  83.3 (28.5)  81.5 (27.5)  --
%    50% work    88.7 (31.5)  85.1 (29.5)  83.3 (28.5)  81.5 (27.5)
%    25% work    90.5 (32.5)  87.8 (31.0)  86.0 (30.0)  85.1 (29.5)
%
% Not acclimatized (ACGIH Action Limit)
%   Work/rest    Light        Moderate     Heavy        Very Heavy
%   100% work    81.5 (27.5)  77.0 (25.0)  72.5 (22.5)  --
%    75% work    84.2 (29.0)  79.7 (26.5)  76.1 (24.5)  --
%    50% work    86.0 (30.0)  82.4 (28.0)  79.7 (26.5)  77.0 (25.0)
%    25% work    87.8 (31.0)  84.2 (29.0)  82.4 (28.0)  79.7 (26.5)
%
% "--": no limit given (very heavy work is not done for 75-100% of a shift).
%
% Inputs
%
% WBGT = Wet Bulb Globe Temperature (degrees F by default; degrees C with
%        'Units','C').  Any size array.  Note that calcwbgt.m returns
%        degrees C, so either pass 'Units','C' or convert (WBGT*9/5+32).
%
% acclimatized = 0:  not acclimatized
%              = 1:  acclimatized
%
% worktype     = 1:  Light
%              = 2:  Moderate
%              = 3:  Heavy
%              = 4:  Very Heavy
%
% adder (optional; an adder, based on clothing type, added to WBGT)
%   Default value is 0.  In degrees F by default; with 'Units','C' the
%   adder is in degrees C (so multiply the values below by 5/9).
%   A scalar, or an array the same size as WBGT.
%              = 0:  Work clothes (long sleeve shirt and pants)
%              = 6:  Cloth (woven material) overalls
%              = 9:  Double-cloth overalls
%              = 1:  SMS polypropylene coveralls
%              = 2:  Polyolefin coveralls
%   These are as in the original (2017) version.  They appear to be older
%   ACGIH clothing adjustment values (about 3.5, 5, 0.5 and 1 degrees C)
%   converted to degrees F and rounded.  They have not been re-verified:
%   check them against Davis et al. (2016) and the current ACGIH Clothing
%   Adjustment Values (CAVs), which differ for some ensembles.
%
% Optional name-value input
% 'Units' = 'F' (default) or 'C': units of WBGT and adder.
%           With 'C', the inputs are converted to degrees F
%           (WBGT*9/5+32 and adder*9/5) and compared with the table above.
%
% Outputs (same size as WBGT)
% workdemand = percentage of work schedule vs rest schedule
% = 100:  Workers can work safely 100% of their shift
% =  75:  Workers need to rest for 25% of their shift
% =  50:  Workers need to rest for 50% of their shift
% =  25:  Workers need to rest for 75% of their shift
% =   0:  Workers cannot work
% = NaN:  WBGT (or adder) is NaN
%
% Examples
%   heatstress(83,0,2)                % 25 (83 F, non-acclimatized, moderate)
%   heatstress(28.3,0,2,'Units','C')  % 25 (the same in degrees C)
%   heatstress([75 80 85 90],1,1,6)   % [100 75 0 0] (cloth overalls, +6 F)
%   heatstress(calcwbgt(Ta,RH,S,p,V,tau,zday,lat),1,3,'Units','C')
%
% References
% Davis, D. W., et al. (2016), Hanford Climate Vulnerability Assessment and
%   Resilience Action Plan, MSA-1604467, Mission Support Alliance,
%   Richland, WA.
% American Conference of Governmental Industrial Hygienists (ACGIH),
%   TLVs and BEIs: Threshold Limit Values for Chemical Substances and
%   Physical Agents and Biological Exposure Indices, "Heat Stress and
%   Strain" (screening criteria for TLV and Action Limit), ACGIH,
%   Cincinnati, OH (annual).

arguments
    WBGT {mustBeNumeric}
    acclimatized (1,1) {mustBeNumericOrLogical,mustBeMember(acclimatized,[0 1])}
    worktype (1,1) {mustBeNumeric,mustBeMember(worktype,[1 2 3 4])}
    adder {mustBeNumeric} = 0
    opts.Units {mustBeTextScalar} = 'F'
end

units=upper(char(opts.Units));
if ~any(strcmp(units,{'F','C'}))
    error('heatstress:units','Units must be ''F'' or ''C''.');
end

% Adding adder to WBGT (or default value of 0), in degrees F
if strcmp(units,'C')
    T=(WBGT*9/5+32)+adder*9/5;
else
    T=WBGT+adder;
end

accwork=[85.1 81.5 78.8 NaN ; 86.9 83.3 81.5 NaN ; 88.7 85.1 83.3 81.5 ; 90.5 87.8 86.0 85.1];
nonaccwork=[81.5 77.0 72.5 NaN ; 84.2 79.7 76.1 NaN ; 86.0 82.4 79.7 77.0 ; 87.8 84.2 82.4 79.7];

% is the worker acclimatized or not?
if acclimatized==0
    thework=nonaccwork;
else
    thework=accwork;
end

% isolating the type of work
thework=thework(:,worktype);

demands=[100 75 50 25];

% The result is the demand of the FIRST row whose limit is >= T (0 if
% none).  Looping from the last row to the first, later (earlier-row)
% assignments overwrite, which gives the first such row.  NaN limits are
% never satisfied.
workdemand=zeros(size(T));
for r=numel(demands):-1:1
    workdemand(T<=thework(r))=demands(r);
end
workdemand(isnan(T))=NaN;
