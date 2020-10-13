function varargout=frostpoint(flag,varargin)

% [output1,[output2]]=frostpoint(flag,input1,input2)
%
% Original code by Susan Solomon
% Modifications by Ben Kravitz (ben.kravitz.work@gmail.com) 16 February 2012
%
% This program can be used for two purposes:
% flag=1:  Calculate the frost point (output1) at a given pressure (input1)
%          for given mass mixing ratios (input2) of water.
% flag=2:  Calculate the saturation mixing ratio over water (output1) and
%          over ice (output2) for a given pressure (input1) and
%          temperature (input2).
% Input units:  p (mb), T (K), water (kg/kg)
%
% References
%
% Flatau, P. J., R. L. Walko, and W. R. Cotton (1992), Polynomial fits to
%     saturation vapor pressure, J. Appl. Met., 31, 1507-1513.
% Wexler, A. (1976), Vapor pressure formulation for water in range 0 to 100°C, A
%    revision, J. Res. Natl. Bur. Stand., 80A, 775-785.
% Wexler, A. (1977), Vapor pressure formulation for ice, J. Res. Natl. Bur.
%    Stand., 81A, 5-20.

if flag==1;
    p=varargin{1};
    water=varargin{2};
    frost=(1.1313*log(p*water)+6182)./(24.45-log(p*water));
    varargout={frost};
elseif flag==2
    p=varargin{1};
    T=varargin{2};
    % Clausius Clapeyron - log of saturation vapor pressure over water (Pa)
    lnewsat=log(611.73)+(2.453e6/461.51)*(1./273.16-1./T);
    % Saturation mixing ratio (kg/kg)
    wsw=0.622*exp(lnewsat)./(p-exp(lnewsat));

    % Saturation Mixing Ratio Over Ice (kg/kg)
    % Approximation from Flatau et al. (1992)
    % Based on polynomial fits to Wexler (1976, 1977)
    wsi=6.11176750+0.443986062*(T-273.15)+...
        0.143053301e-1*(T-273.15).^2+...
        0.265027242e-3*(T-273.15).^3+...
        0.302246994e-5*(T-273.15).^4+...
        0.203886313e-7*(T-273.15).^5+...
        0.638780966e-10*(T-273.15).^6;
    
    varargout(1)={wsw};
    varargout(2)={wsi};
end

% If you want the values to write out to a file, feel free to
% uncomment these lines.
%
% fid=fopen('frost.out','w');
% fprintf(fid,['Frost point for ' num2str(p) ' mb\n']);
% fprintf(fid,'   Water    Frost Point\n');
% for k=1:length(water);
%     fprintf(fid,'%1.4e   %3.4f\n',water(k),frost(k));
% end
% fprintf(fid,'\n\n');
% fprintf(fid,'Saturation Mixing Ratio\n');
% fprintf(fid,'Temperature  |  Over H2O  |  Over Ice\n');
% for k=1:length(T);
%     fprintf(fid,'  %3.4f     %1.4e   %1.4e\n',T(k),wsw(k),wsi(k));
% end
% fclose(fid);
