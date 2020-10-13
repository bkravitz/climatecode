function T=equilibrium(L,epsilon,S0,alpha,accuracy)

% function T=equilibrium(L,epsilon,S0,alpha,accuracy)
% A simple radiative equilibrium energy balance model with an
% arbitrary number of layers (at least two).
%
% Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Originally written 12 October 2011
% Functionalized 13 February 2018
% Fixed math error 5 April 2020
%
% Inputs:
%   L = number of layers (must be at least 2:  surface and TOA) [default = 15]
%   epsilon = array of emissivities for each layer (1=opaque, 0=transparent) [default = 0.7*ones(1,L)]
%   S0 = solar constant (W/m2) [default=1366]
%   alpha = surface albedo [default=0.3]
%   accuracy = Defines when the model quits iterating [default=0.01]
%
% Output:  T = temperature at each layer (Kelvin)

%%%%%%%%%%%%%% model begins here %%%%%%%%%%%%%%
%%% modify anything below at your own peril %%%

%% constants
sigma=5.67e-8; % Stefan-Boltzmann constant
T=zeros(1,L); % initializing temperatures array
surf_sw=(S0/4)*(1-alpha);
if T(1)==0;
    T(1)=(surf_sw/sigma)^(1/4); % initial surface temperature (shortwave only)
    % the model is relatively insensitive to the choice of initial surface temperature
    % so I just picked something reasonable
end
T2=zeros(1,L); % initializing comparison array

%% iterate over everything and stop when temperature change is small
while norm(T-T2)>accuracy;
    T2=T;
    
    rads=epsilon.*sigma.*T.^4;

    %% surface (layer 1)
    surf_sw=(S0/4)*(1-alpha);

    terms1=zeros(1,L); % contributions of each layer to layer 1
    terms1(2:L)=rads(2:L);
    for n=2:L;
        counter=n;
        while counter>2;
            terms1(n)=terms1(n)*(1-epsilon(counter-1));
            counter=counter-1;
        end
    end
    %surf_lw=sum(terms1)*epsilon(1);
    surf_lw=sum(terms1);
    T(1)=((surf_sw+surf_lw)/sigma)^(1/4); % adjusting surface temperature to include longwave

    %% layer n (2 to L-1) if these layers exist (i.e., if L>2)
    if L>2;
        for n=2:L-1;
            termsn=zeros(1,L); % contributions of each layer to layer n
            termsn(n+1:L)=rads(n+1:L); % layers above
            for k=n+1:L;
                counter=k;
                while counter>n+1;
                    termsn(k)=termsn(k)*(1-epsilon(counter-1));
                    counter=counter-1;
                end
            end
            termsn(1:n-1)=rads(1:n-1); % layers below
            for k=1:n-1;
                counter=k;
                while counter<n-1;
                   termsn(k)=termsn(k)*(1-epsilon(counter+1));
                   counter=counter+1;
                end
            end
            %layer_n_lw=sum(termsn)*epsilon(n);
            layer_n_lw=sum(termsn);
            T(n)=(0.5*layer_n_lw/sigma)^(1/4);
        end
    end

    %% layer L (TOA)
    termsl=zeros(1,L); % contributions of each layer to layer L
    termsl(1:L-1)=rads(1:L-1);
    for k=1:L-1;
        counter=k;
        while counter<L-1;
            termsl(k)=termsl(k)*(1-epsilon(counter+1));
            counter=counter+1;
        end
    end
    %layer_L_lw=sum(termsl)*epsilon(L);
    layer_L_lw=sum(termsl);
    T(L)=(0.5*layer_L_lw/sigma)^(1/4);
end
