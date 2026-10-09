function T=equilibrium(L,epsilon,S0,alpha,accuracy)

% function T=equilibrium(L,epsilon,S0,alpha,accuracy)
% A simple radiative equilibrium energy balance model with an
% arbitrary number of layers (at least two).  Layer 1 is the surface,
% which absorbs all shortwave and longwave reaching it and emits as a
% blackbody.  Layers 2 to L are gray longwave absorbers/emitters that are
% transparent to shortwave and emit eps*sigma*T^4 up and down.
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Originally written 12 October 2011
% Functionalized 13 February 2018
% Fixed math error 5 April 2020
% Surface made a blackbody emitter (energy conservation) 8 October 2026
%
% Inputs:
%   L = number of layers (must be at least 2:  surface and TOA) [default = 15]
%   epsilon = array of emissivities for each layer (1=opaque, 0=transparent) [default = 0.7*ones(1,L)]
%             (scalar is expanded to all layers; epsilon(1) is ignored
%             because the surface is a blackbody)
%   S0 = solar constant (W/m2) [default=1366]
%   alpha = surface albedo [default=0.3]
%   accuracy = Defines when the model quits iterating:  stop when the
%              norm of the change in T between iterations (K) is below this [default=0.01]
%
% Output:  T = 1 x L array of temperature at each layer (Kelvin), T(1) = surface

%%%%%%%%%%%%%% model begins here %%%%%%%%%%%%%%
%%% modify anything below at your own peril %%%

%% defaults and input checks
if nargin<1 || isempty(L), L=15; end
if nargin<2 || isempty(epsilon), epsilon=0.7*ones(1,L); end
if nargin<3 || isempty(S0), S0=1366; end
if nargin<4 || isempty(alpha), alpha=0.3; end
if nargin<5 || isempty(accuracy), accuracy=0.01; end
if ~isscalar(L) || L<2 || L~=round(L)
    error('equilibrium:badL','L must be an integer >= 2.');
end
if isscalar(epsilon), epsilon=epsilon*ones(1,L); end
if numel(epsilon)~=L
    error('equilibrium:badEpsilon','epsilon must have L elements.');
end
epsilon=reshape(epsilon,1,L); % row vector, to match T

%% constants
sigma=5.67e-8; % Stefan-Boltzmann constant
T=zeros(1,L); % initializing temperatures array
surf_sw=(S0/4)*(1-alpha);
if T(1)==0
    T(1)=(surf_sw/sigma)^(1/4); % initial surface temperature (shortwave only)
    % the model is relatively insensitive to the choice of initial surface temperature
    % so I just picked something reasonable
end
T2=zeros(1,L); % initializing comparison array

%% iterate over everything and stop when temperature change is small
while norm(T-T2)>accuracy
    T2=T;
    
    rads=epsilon.*sigma.*T.^4; % longwave emitted by each layer (each direction)
    rads(1)=sigma*T(1)^4; % surface emits as a blackbody, consistent with T(1) below

    %% surface (layer 1)
    surf_sw=(S0/4)*(1-alpha);

    terms1=zeros(1,L); % contributions of each layer to layer 1 (attenuated by layers in between)
    terms1(2:L)=rads(2:L);
    for n=2:L
        counter=n;
        while counter>2
            terms1(n)=terms1(n)*(1-epsilon(counter-1));
            counter=counter-1;
        end
    end
    %surf_lw=sum(terms1)*epsilon(1);
    surf_lw=sum(terms1);
    T(1)=((surf_sw+surf_lw)/sigma)^(1/4); % surface absorbs all sw+lw and emits sigma*T^4

    %% layer n (2 to L-1) if these layers exist (i.e., if L>2)
    if L>2
        for n=2:L-1
            termsn=zeros(1,L); % contributions of each layer to layer n
            termsn(n+1:L)=rads(n+1:L); % layers above
            for k=n+1:L
                counter=k;
                while counter>n+1
                    termsn(k)=termsn(k)*(1-epsilon(counter-1));
                    counter=counter-1;
                end
            end
            termsn(1:n-1)=rads(1:n-1); % layers below
            for k=1:n-1
                counter=k;
                while counter<n-1
                   termsn(k)=termsn(k)*(1-epsilon(counter+1));
                   counter=counter+1;
                end
            end
            %layer_n_lw=sum(termsn)*epsilon(n);
            layer_n_lw=sum(termsn);
            T(n)=(0.5*layer_n_lw/sigma)^(1/4); % absorbs eps*lw, emits 2*eps*sigma*T^4
        end
    end

    %% layer L (TOA)
    termsl=zeros(1,L); % contributions of each layer to layer L
    termsl(1:L-1)=rads(1:L-1);
    for k=1:L-1
        counter=k;
        while counter<L-1
            termsl(k)=termsl(k)*(1-epsilon(counter+1));
            counter=counter+1;
        end
    end
    %layer_L_lw=sum(termsl)*epsilon(L);
    layer_L_lw=sum(termsl);
    T(L)=(0.5*layer_L_lw/sigma)^(1/4);
end
