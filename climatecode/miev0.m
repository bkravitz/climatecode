function [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,pmom,spike]=...
         miev0(xx,crefin,perfct,mimcut,anyang,numang,xmu,nmom,...
         ipolzn,momdim)

% function [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,pmom,spike]=
%          miev0(xx,crefin,perfct,mimcut,anyang,numang,xmu,nmom,...
%          ipolzn,momdim)
%
%    MIEV0
%    Computes Mie scattering and extinction efficiencies; asymmetry
%    factor;  forward- and backscatter amplitude;  scattering
%    amplitudes vs. scattering angle for incident polarization parallel
%    and perpendicular to the plane of scattering;
%    coefficients in the Legendre polynomial expansions of either the
%    unpolarized phase function or the polarized phase matrix;
%    some quantities needed in polarized radiative transfer;  and
%    information about whether or not a resonance has been hit.
%
%    MATLAB translation of Warren Wiscombe's Fortran MIEV0 (Matlabized by
%    Ben Kravitz, May 2012; subroutines are in this folder).
%    References:
%      Wiscombe, W.J., 1979: Mie Scattering Calculations: Advances in
%        Technique and Fast, Vector-Speed Computer Codes. NCAR Tech. Note
%        NCAR/TN-140+STR (edited/revised 1996).
%      Wiscombe, W.J., 1980: Improved Mie scattering algorithms. Appl. Opt.,
%        19, 1505-1509.
%    Single sphere, single size parameter.  All efficiencies are
%    dimensionless.  For very small particles (xx*max(1,|m|) <= 0.1, or
%    xx <= 0.1 if perfct==1) the 2-term small-particle limits small1/small2
%    are used instead of the full series.
%
%    Bug fixes 8 October 2026 (with Claude): ipolzn digits are now parsed
%      directly (any order/combination; invalid digits give an error)
%      instead of via a lookup table; xmu length checked; outputs s1, s2,
%      tforw, tback, pmom are defined on every path (pmom=[] if nmom==0).
%      Results-changing fix in lpcoef: S21 (ipolzn digit 3) moments were
%      computed from C*conj(C) instead of C*conj(D).  pmom is now
%      (nmom+1) x (number of phase quantities); momdim is no longer used.
%      spike is now 1 (no resonance) in the small-particle limits instead
%      of 0, as documented below; global variables removed.
%    Changes 8 October 2026 (with Claude): MIEV0 consolidated into this
%      single file.  The former separate files small1, small2, biga,
%      confra, lpcoef, lpco1t, lpco2t and sq are now local functions below
%      (unchanged); the statement functions f1 and f2 are defined inside
%      biga, as in the Fortran; f3 (never used) was removed.  Results are
%      identical to the multi-file version.
%
%    Inputs:
%      xx     - mie size parameter (2*pi*radius/wavelength)
%      crefin - complex refractive index (imaginary part can be
%               positive or negative)
%      perfct - If 1, assumes refractive index is infinite.
%               (Common case:  perfct=0)
%      mimcut - Positive value below which imaginary refractive index
%               is regarded as 0 (speeds up computation).
%      anyang - If 1, any angles may be input through xmu.  If 0,
%               angles are monotone increasing and symmetric about 90
%               degrees.
%      numang - Number of angles at which scattering amplitudes are to
%               be calculated (0 = none).
%      xmu    - Vector (at least numang elements) of cosines of angles at which s1 and s2 are to be
%               evaluated.  If anyang==0, then the angles must be
%               monotone increasing and symmetric about 90 degrees.
%               If numang is odd and anyang==0, 90 degrees must be
%               among the angles.
%      nmom   - Highest Legendre moment of PMOM to calculate.  Numbering
%               starts from 0.
%      ipolzn - POSITIVE:  Compute Legendre moments for the Mueller matrix
%                   elements determined by the digits of ipolzn:
%                   1 --> M1, 2 --> M2, 3 --> S21, 4 --> D21.  Sample input
%                   ipolzn=14.
%               0:  Compute Legendre moments pmom for the unpolarized,
%                   unnormalized phase function.
%               NEGATIVE:  Computes Legendre moments for the Sekera phase
%                   quantities determined by the digits of abs(ipolzn):
%                   1 --> R1, 2 --> R2, 3 --> R3, 4 --> R4.  Sample input
%                   ipolzn=-14.
%               Not used if nmom=0;
%      momdim - Kept for compatibility with the Fortran calling sequence;
%               no longer used (pmom is sized from nmom and ipolzn).
%
%    Outputs:
%      qext   - Extinction efficiency factor
%      qsca   - Scattering efficiency factor
%      gqsc   - Asymmetry factor times scattering efficiency
%      s1(n)  - Complex Mie scattering amplitudes at angles xmu(n)
%      s2(n)    (1 x numang; s1 perpendicular, s2 parallel polarization)
%      sforw  - Complex forward-scattering amplitude s1 at 0 degrees
%      sback  - Complex backscattering amplitude s1 at 180 degrees
%      tforw  - 1x2 complex [T1 T2] at 0 degrees, needed for polarized
%               radiative transfer (see Wiscombe 1979)
%      tback  - 1x2 complex [T1 T2] at 180 degrees
%      spike  - Smallest Mie-coefficient denominator over terms n > xx.
%               Values below about 0.3 indicate a resonance (ripple
%               spike); 1 when no spike is detected, including the
%               small-particle limits (no resonances there).
%      pmom(m+1,np) - moments m=0 to nmom of unnormalized np-th phase
%               quantity; (nmom+1) x (number of phase quantities),
%               row 1 = moment 0.  Empty if nmom==0.
%
%    Local functions (below):  small1, small2, biga, confra, lpcoef,
%      lpco1t, lpco2t, sq

if numel(xmu)<numang
    error('miev0:numang','xmu must have at least numang elements.');
end
% calcmo(k)==1 if phase quantity k is requested; npquan = highest k
npquan=0;
calcmo=zeros(1,4);
pmom=[];
if nmom~=0
    if ipolzn~=0
        % digits of abs(ipolzn) select the phase quantities (1-4)
        ipdig=num2str(abs(ipolzn))-'0';
        if any(ipdig<1 | ipdig>4)
            error('miev0:ipolzn','Each digit of ipolzn must be 1, 2, 3, or 4.');
        end
        calcmo(ipdig)=1;
        npquan=max(ipdig);
    end
end

s1=zeros(1,numang);
s2=zeros(1,numang);
countt=0;  % set to 1 once a small-particle limit has supplied the answer
if (perfct==1) && (xx<=0.1)
    [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,lita,litb]=small1(xx,numang,xmu);
    ntrm=2;
    countt=1;
    spike=1; % no resonance spikes in the small-particle limit
end
noabs=1;
if countt==0
    if perfct==0
        cior=crefin;
        % internally use the Im(m)<=0 convention; outputs are conjugated
        % back at the end if the user gave Im(m)>0
        if imag(cior)>0, cior=conj(cior); end
        mre=real(cior);
        mim=-imag(cior);
        noabs=double(mim<=mimcut);
        cioriv=1./cior;
        rioriv=1./mre;
        if xx*max([1 abs(cior)])<=0.1
            templog=double(mim>mimcut);
            [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,lita,litb]=small2(xx,cior,templog,numang,xmu);
            ntrm=2;
            countt=1;
            spike=1; % no resonance spikes in the small-particle limit
        end
    end
end
if countt==0
    nangd2=floor((numang+1)/2);
    yesang=double(numang>0);
    % number of terms in the Mie series (Wiscombe 1980 criterion)
    if xx<=8
        ntrm=floor(xx+4*xx.^(1/3)+1);
    elseif xx<4200
        ntrm=floor(xx+4.05*xx.^(1/3)+2);
    else
        ntrm=floor(xx+4*xx.^(1/3)+2);
    end
    lita=zeros(1,ntrm);
    litb=zeros(1,ntrm);
    if perfct==0
        bigaout=biga(cior,xx,ntrm,noabs,yesang);
        if noabs==1
            rbiga=bigaout;
        else
            cbiga=bigaout;
        end
    end
    xinv=1./xx;
    psinm1=sin(xx);
    chinm1=cos(xx);
    psin=psinm1*xinv-chinm1;
    chin=chinm1*xinv+psinm1;
    zetnm1=psinm1+1i*chinm1;
    zetn=psin+1i*chin;
    anm1=0;
    bnm1=0;
    if anyang==1
        pinm1=zeros(1,numang);
        pin=ones(1,numang);
        sp=zeros(1,numang);
        sm=zeros(1,numang);
    else
        pinm1=zeros(1,nangd2);
        pin=ones(1,nangd2);
        sp=zeros(1,nangd2);
        sm=zeros(1,nangd2);
        sps=zeros(1,nangd2);
        sms=zeros(1,nangd2);
    end

    qsca=0;
    gqsc=0;
    sforw=0;
    sback=0;
    tforw=zeros(1,2);
    tback=zeros(1,2);
    csum1=0;
    csum2=0;

    % loop to sum mie series
    mm=1;
    spike=1;
    for n=1:ntrm
        fn=n;
        rn=1/fn;
        np1dn=1+rn;
        twonp1=2*n+1;
        coeff=twonp1/(fn*(n+1));
        tcoef=twonp1*(fn*(n+1));
        if perfct==1 % totally reflecting case
            an=((fn*xinv)*psin-psinm1)/((fn*xinv)*zetn-zetnm1);
            bn=psin/zetn;
        elseif noabs==1 % no absorption case
            cdenan=(rioriv*rbiga(n)+(fn*xinv))*zetn-zetnm1;
            an=((rioriv*rbiga(n)+(fn*xinv))*psin-psinm1)/cdenan;
            cdenbn=(mre*rbiga(n)+(fn*xinv))*zetn-zetnm1;
            bn=((mre*rbiga(n)+(fn*xinv))*psin-psinm1)/cdenbn;
        else
            cdenan=(cioriv*cbiga(n)+(fn*xinv))*zetn-zetnm1;
            cdenbn=(cior*cbiga(n)+(fn*xinv))*zetn-zetnm1;
            an=((cioriv*cbiga(n)+(fn*xinv))*psin-psinm1)/cdenan;
            bn=((cior*cbiga(n)+(fn*xinv))*psin-psinm1)/cdenbn;
            qsca=qsca+twonp1*(sq(an)+sq(bn));
        end
        lita(n)=an;
        litb(n)=bn;
        % resonance (spike) detection: very unequal, small denominators
        if (perfct==0) && (n>xx)
            denan=abs(cdenan);
            denbn=abs(cdenbn);
            ratio=denan/denbn;
            if (ratio<=0.2) || (ratio>=5.0)
                spike=min([spike denan denbn]);
            end
        end
        sforw=sforw+twonp1*(an+bn);
        csum1=csum1+tcoef*(an-bn);
        sback=sback+(mm*twonp1)*(an-bn);
        csum2=csum2+(mm*tcoef)*(an+bn);
        gqsc=gqsc+(fn-rn)*real(anm1*conj(an)+bnm1*conj(bn))+coeff*real(an*conj(bn));
        if yesang==1
            anp=coeff*(an+bn);
            bnp=coeff*(an-bn);
            if anyang==1
                for j=1:numang
                    rtmp=(xmu(j)*pin(j))-pinm1(j);
                    taun=fn*rtmp-pinm1(j);
                    sp(j)=sp(j)+anp*(pin(j)+taun);
                    sm(j)=sm(j)+bnp*(pin(j)-taun);
                    pinm1(j)=pin(j);
                    pin(j)=(xmu(j)*pin(j))+np1dn*rtmp;
                end
            else
                anpm=mm*anp;
                bnpm=mm*bnp;
                for j=1:nangd2
                    rtmp=(xmu(j)*pin(j))-pinm1(j);
                    taun=fn*rtmp-pinm1(j);
                    sp(j)=sp(j)+anp*(pin(j)+taun);
                    sms(j)=sms(j)+bnpm*(pin(j)+taun);
                    sm(j)=sm(j)+bnp*(pin(j)-taun);
                    sps(j)=sps(j)+anpm*(pin(j)-taun);
                    pinm1(j)=pin(j);
                    pin(j)=(xmu(j)*pin(j))+np1dn*rtmp;
                end
            end
        end
        mm=-mm;
        anm1=an;
        bnm1=bn;
        zet=(twonp1*xinv)*zetn-zetnm1;
        zetnm1=zetn;
        zetn=zet;
        psinm1=psin;
        psin=real(zetn);
    end

    qext=2./(xx.^2)*real(sforw);
    if (perfct==1) || (noabs==1)
        qsca=qext;
    else
        qsca=2./xx.^2*qsca;
    end
    gqsc=4./xx.^2*gqsc;
    
    sforw=0.5*sforw;
    sback=0.5*sback;
    tforw(1)=0.5*sforw-0.125*csum1;
    tforw(2)=0.5*sforw+0.125*csum1;
    tback(1)=-0.5*sback+0.125*csum2;
    tback(2)=0.5*sback+0.125*csum2;
    if yesang==1
        if anyang==1
            s1(1:numang)=0.5*(sp(1:numang)+sm(1:numang));
            s2(1:numang)=0.5*(sp(1:numang)-sm(1:numang));
        else
            s1(1:nangd2)=0.5*(sp(1:nangd2)+sm(1:nangd2));
            s2(1:nangd2)=0.5*(sp(1:nangd2)-sm(1:nangd2));
            for kka=1:nangd2
                s1(numang+1-kka)=0.5*(sps(kka)+sms(kka));
                s2(numang+1-kka)=0.5*(sps(kka)-sms(kka));
            end
        end
    end
end

if nmom>0
    pmom=lpcoef(ntrm,nmom,ipolzn,momdim,calcmo,npquan,lita,litb);
end
if imag(crefin)>0
    sforw=conj(sforw);
    sback=conj(sback);
    tforw(1:2)=conj(tforw(1:2));
    tback(1:2)=conj(tback(1:2));
    s1(1:numang)=conj(s1(1:numang));
    s2(1:numang)=conj(s2(1:numang));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,a,b]=small1(xx,numang,xmu)

% function [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,a,b]=small1(xx,numang,xmu)
% Small-particle limit of Mie quantities in totally reflecting
%  limit (Mie series truncated after 2 terms)
% Subroutine used in MIEV0 (Fortran routine SMALL1); used when perfct==1
%  and xx <= 0.1
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 1 May 2012
% Code cleanup 8 October 2026 (with Claude): a, b, tforw, tback
%  preallocated; results unchanged.
%
% Inputs:  xx, numang, xmu - as in miev0
% Output:  qext, qsca, gqsc, sforw, sback, s1, s2, tforw, tback - as in
%          miev0; a, b - first two Mie coefficients a_n, b_n (1x2 complex)
% Functions used:  sq

% preallocations
s1=zeros(1,numang);
s2=zeros(1,numang);
a=zeros(1,2);
b=zeros(1,2);
tforw=zeros(1,2);
tback=zeros(1,2);

a(1)=(0+1i*(2/3)*(1-0.2*xx.^2))./(1-0.5*xx.^2+1i*(2/3)*xx.^3);
b(1)=(0-1i*(1-0.1*xx.^2)./3)./(1+0.5*xx.^2-1i*xx.^3/3);
a(2)=0+1i*xx.^2/30;
b(2)=0-1i*xx.^2/45;
qsca=6*xx.^4*(sq(a(1))+sq(b(1))+(5/3)*(sq(a(2))+sq(b(2))));
qext=qsca;
gqsc=6*xx.^4*real(a(1)*conj(a(2)+b(1))+(b(1)+(5/9)*a(2))*conj(b(2)));
rtmp=1.5*xx^3;
sforw=rtmp*(a(1)+b(1)+(5/3)*(a(2)+b(2)));
sback=rtmp*(a(1)-b(1)-(5/3)*(a(2)-b(2)));
tforw(1)=rtmp*(b(1)+(5/3)*(2*b(2)-a(2)));
tforw(2)=rtmp*(a(1)+(5/3)*(2*a(2)-b(2)));
tback(1)=rtmp*(b(1)-(5/3)*(2*b(2)+a(2)));
tback(2)=rtmp*(a(1)-(5/3)*(2*a(2)+b(2)));

s1(1:numang)=rtmp.*(a(1)+b(1).*xmu(1:numang)+(5/3).*(a(2).*xmu(1:numang)+b(2).*(2.*xmu(1:numang).^2-1)));
s2(1:numang)=rtmp.*(b(1)+a(1).*xmu(1:numang)+(5/3).*(b(2).*xmu(1:numang)+a(2).*(2.*xmu(1:numang).^2-1)));

% Recover actual Mie coefficients
a(1)=xx.^3*a(1);
a(2)=xx.^3*a(2);
b(1)=xx.^3*b(1);
b(2)=xx.^3*b(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,a,b]=small2(xx,cior,calcqe,numang,xmu)

% function [qext,qsca,gqsc,sforw,sback,s1,s2,tforw,tback,a,b]=small2(xx,cior,calcqe,numang,xmu)
% Small-particle limit of Mie quantities for general refractive index
%  (Mie series truncated after 2 terms)
% Subroutine used in MIEV0 (Fortran routine SMALL2); used when
%  xx*max(1,|cior|) <= 0.1
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 2 May 2012
% Code cleanup 8 October 2026 (with Claude): a, b, tforw, tback
%  preallocated; results unchanged.
%
% Inputs:  xx, numang, xmu - as in miev0
%          cior   - complex refractive index (Im(cior) <= 0 convention)
%          calcqe - 1 to compute qext separately (absorbing particle);
%                   0 sets qext=qsca
% Output:  qext, qsca, gqsc, sforw, sback, s1, s2, tforw, tback - as in
%          miev0; a, b - first two Mie coefficients a_n, b_n (1x2
%          complex; b(2)=0 in this approximation)
% Functions used:  sq

% preallocations
s1=zeros(1,numang);
s2=zeros(1,numang);
a=zeros(1,2);
b=zeros(1,2);
tforw=zeros(1,2);
tback=zeros(1,2);

ciorsq=cior.^2;
ctmp=(0+1i*(2/3))*(ciorsq-1);
a(1)=ctmp*(1-0.1*xx.^2+(ciorsq/350+(1/280))*xx.^4)./(ciorsq+ ...
    2+(1-0.7*ciorsq)*xx.^2-(ciorsq.^2/175-0.275*ciorsq+0.25)*xx.^4+xx.^3*ctmp*(1-0.1*xx.^2));
b(1)=(xx.^2/30)*ctmp*(1+(ciorsq/35-(1/14))*xx.^2)./(1-(ciorsq/15-(1/6))*xx.^2);
a(2)=(0.1*xx.^2)*ctmp*(1-xx.^2/14)./(2*ciorsq+3-(ciorsq/7-0.5)*xx.^2);
qsca=(6*xx.^4)*(sq(a(1))+sq(b(1))+(5/3)*sq(a(2)));
qext=qsca;
if calcqe==1;qext=6*xx*real(a(1)+b(1)+(5/3)*a(2));end
gqsc=(6*xx.^4)*real(a(1)*conj(a(2)+b(1)));
rtmp=1.5*xx.^3;
sforw=rtmp*(a(1)+b(1)+(5/3)*a(2));
sback=rtmp*(a(1)-b(1)-(5/3)*a(2));
tforw(1)=rtmp*(b(1)-(5/3)*a(2));
tforw(2)=rtmp*(a(1)+2*(5/3)*a(2));
tback(1)=tforw(1);
tback(2)=rtmp*(a(1)-2*(5/3)*a(2));
s1(1:numang)=rtmp.*(a(1)+(b(1)+(5/3).*a(2)).*xmu(1:numang));
s2(1:numang)=rtmp.*(b(1)+a(1).*xmu(1:numang)+(5/3).*a(2).*(2.*xmu(1:numang).^2-1));
a(1)=xx.^3*a(1);
a(2)=xx.^3*a(2);
b(1)=xx.^3*b(1);
b(2)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output=biga(cior,xx,ntrm,noabs,yesang)

% function output=biga(cior,xx,ntrm,noabs,yesang)
% Calculates logarithmic derivatives of J-Bessel function
%  A_n(m*x) = psi_n'(m*x)/psi_n(m*x), n=1..ntrm, by up- or down-recurrence
%  (choice per Wiscombe 1980, Appl. Opt. 19, 1505-1509)
% Subroutine used in MIEV0 (Fortran routine BIGA)
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 1 May 2012
% Code cleanup 8 October 2026 (with Claude): syntax only, results unchanged.
%
% Inputs:  cior   - complex refractive index (Im(cior) <= 0 convention)
%          xx     - size parameter
%          ntrm   - number of terms in the Mie series
%          noabs  - 1 if no absorption (real arithmetic), else 0
%          yesang - 1 if scattering angles are requested
% Output:  rbiga (real, 1 x ntrm) if noabs==1, else cbiga (complex, 1 x ntrm)
% Functions used:  confra

% Statement functions (Wiscombe 1980 fits in terms of the real refractive
% index): if mim*xx is below f1 (no angles) or f2 (angles requested),
% up-recurrence is stable.
f1=@(mre) -8+(mre.^2).*(26.22+mre.*(-0.4474+(mre.^3).*(0.00204-0.000175*mre)));
f2=@(mre) 3.9+mre.*(-10.8+13.78*mre);

% pre-allocations
rbiga=zeros(1,ntrm);
cbiga=zeros(1,ntrm);

% Decide whether BigA can be calculated by up-recurrence
mre=real(cior);
mim=abs(imag(cior));
if (mre<1) || (mre>10) || (mim>10)
    down=1;
else
    if yesang==1
        down=1;
        if (mim*xx<f2(mre))
            down=0;
        end
    else
        down=1;
        if (mim*xx<f1(mre))
            down=0;
        end
    end
end

zinv=1./(cior.*xx);
rezinv=1./(mre.*xx);

if down==1 % Downward recurrence for BigA
    % Compute initial high-order BigA using Lentz method
    ctmp=confra(ntrm,zinv);
    if noabs==1  % no-absorption case
        rbiga(ntrm)=real(ctmp);
        for n=ntrm:-1:2
            rbiga(n-1)=(n*rezinv)-1./((n*rezinv)+rbiga(n));
        end
    else  % absorptive case
        cbiga(ntrm)=ctmp;
        for n=ntrm:-1:2
            cbiga(n-1)=(n*zinv)-1./((n*zinv)+cbiga(n));
        end
    end
else % Upward recurrence for BigA
    if noabs==1  % no-absorption case
        rtmp=sin(mre*xx);
        rbiga(1)=-rezinv+rtmp./(rtmp*rezinv-cos(mre*xx));
        for n=2:ntrm
            rbiga(n)=-(n*rezinv)+1./((n*rezinv)-rbiga(n-1));
        end
    else  % absorptive case
        ctmp=exp(-2*1i*cior*xx);
        cbiga(1)=-zinv+(1-ctmp)./(zinv*(1-ctmp)-1i*(1+ctmp));
        for n=2:ntrm
            cbiga(n)=-(n*zinv)+1./((n*zinv)-cbiga(n-1));
        end
    end
end

if noabs==1
    output=rbiga;
else
    output=cbiga;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output=confra(n,zinv)

% function output=confra(n,zinv)
% Compute Bessel function ratio A-sub-N from its continued fraction using Lentz method
%  (Lentz 1976, Appl. Opt. 15, 668-671); starting value for the
%  down-recurrence in biga
% Subroutine used in MIEV0 (Fortran routine CONFRA)
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 1 May 2012
% Bug fixes 8 October 2026 (with Claude): added an iteration limit
%  (10000) so non-convergence gives an error instead of an infinite loop.
%
% Inputs:  n    - order of A_n
%          zinv - reciprocal of the (complex) argument, 1/(m*x)
% Output:  output - A_n(1/zinv)

eps1=0.01;
eps2=1e-8;

output=(n+1)*zinv;
mm=-1;
kk=2*n+3;
cak=(mm*kk)*zinv;
cdenom=cak;
cnumer=cdenom+1./output;

maxit=10000;  % iteration cap (kount counts terms, as in the Fortran)
kount=1;
convflag=0;
while convflag<1
    kount=kount+1;
    if kount>maxit
        error('confra:noconv','Iteration failed to converge');
    end
    mm=-mm;
    kk=kk+2;
    cak=(mm*kk)*zinv;
    if (abs(cnumer/cak)<=eps1) || (abs(cdenom/cak)<=eps1) % ill-conditioned case
        cntn=cak*cnumer+1;
        cdtd=cak*cdenom+1;
        output=(cntn/cdtd)*output;
        mm=-mm;
        kk=kk+2;
        cak=(mm*kk)*zinv;
        cnumer=cak+cnumer/cntn;
        cdenom=cak+cdenom/cdtd;
        kount=kount+1;
    else
        capt=cnumer/cdenom;
        output=capt*output;
        if (abs(real(capt)-1)>=eps2) || (abs(imag(capt))>=eps2) % check for convergence
            cnumer=cak+1./cnumer;
            cdenom=cak+1./cdenom;
        else
            convflag=1; % values have converged
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmom=lpcoef(ntrm,nmom,ipolzn,momdim,calcmo,npquan,a,b)

% function pmom=lpcoef(ntrm,nmom,ipolzn,momdim,calcmo,npquan,a,b)
% Calculate Legendre polynomial expansion coefficients (also called moments)
%  for phase quantities
% Subroutine used in MIEV0 (Fortran routine LPCOEF).  Uses the double-sum
%  formulas of Dave (1970, Appl. Opt. 9, 1888-1896); see Wiscombe (1979),
%  NCAR/TN-140+STR.
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 3 May 2012
% Bug fixes 8 October 2026 (with Claude): S21 moments (calcmo(3)) now use
%  C*conj(D) as in the Fortran (previously C*conj(C), which gave wrong
%  S21 moments); internal array size maxtrm grows with ntrm so large size
%  parameters no longer overrun the arrays; local "sum" renamed tsum.
%
% Inputs:  ntrm   - number of terms in the Mie series
%          nmom, ipolzn, calcmo, npquan - as set in miev0
%          momdim - unused (kept for the Fortran calling sequence)
%          a, b   - Mie coefficients a_n, b_n (1 x ntrm complex)
% Output:  pmom   - (nmom+1) x max(1,npquan) Legendre moments, row 1 =
%                   moment 0 (unnormalized; see miev0)
% Functions used:  lpco1t, lpco2t

maxtrm=max(1102,ntrm+2);
maxmom=2*maxtrm;
mxmom2=floor(maxmom/2);
maxrcp=4*maxtrm+2;

recip=1./(1:maxrcp);

% Preallocation
pmom=zeros(nmom+1,max(1,npquan));
am=zeros(1,maxtrm+1);
bi=zeros(1,mxmom2+1);
bidel=zeros(1,mxmom2+1);
c=zeros(1,maxtrm);
cm=zeros(1,maxtrm);
cs=zeros(1,maxtrm);
d=zeros(1,maxtrm);
dm=zeros(1,maxtrm);
ds=zeros(1,maxtrm);

if ntrm==1
    pmom=lpco1t(nmom,ipolzn,momdim,calcmo,a,b);
elseif ntrm==2
    pmom=lpco2t(nmom,ipolzn,momdim,calcmo,a,b);
else
    % Calculate Mueller C, D arrays
    cm(ntrm+2)=0;
    dm(ntrm+2)=0;
    cm(ntrm+1)=(1-recip(ntrm+1))*b(ntrm);
    dm(ntrm+1)=(1-recip(ntrm+1))*a(ntrm);
    cm(ntrm)=(recip(ntrm)+recip(ntrm+1))*a(ntrm)+(1-recip(ntrm))*b(ntrm-1);
    dm(ntrm)=(recip(ntrm)+recip(ntrm+1))*b(ntrm)+(1-recip(ntrm))*a(ntrm-1);
    for k=ntrm-1:-1:2
        cm(k)=cm(k+2)-(1+recip(k+1))*b(k+1)+(recip(k)+recip(k+1))*a(k)+(1-recip(k))*b(k-1);
        dm(k)=dm(k+2)-(1+recip(k+1))*a(k+1)+(recip(k)+recip(k+1))*b(k)+(1-recip(k))*a(k-1);
    end
    cm(1)=cm(3)+1.5*(a(1)-b(2));
    dm(1)=dm(3)+1.5*(b(1)-a(2));
    if ipolzn>=0
        c(1:ntrm+2)=(2*(1:ntrm+2)-1).*cm(1:ntrm+2);
        d(1:ntrm+2)=(2*(1:ntrm+2)-1).*dm(1:ntrm+2);
    else
        % Compute Sekera C and D arrays
        cs(ntrm+2)=0;
        ds(ntrm+2)=0;
        cs(ntrm+1)=0;
        ds(ntrm+1)=0;
        for k=ntrm:-1:1
            cs(k)=cs(k+2)+(2*k+1)*(cm(k+1)-b(k));
            ds(k)=ds(k+2)+(2*k+1)*(dm(k+1)-a(k));
        end
        c(1:ntrm+2)=(2*(1:ntrm+2)-1).*cs(1:ntrm+2);
        d(1:ntrm+2)=(2*(1:ntrm+2)-1).*ds(1:ntrm+2);
    end
    if ipolzn<0
        nummom=min([nmom 2*ntrm-2]);
    else
        nummom=min([nmom 2*ntrm]);
    end
    % Loop over moments
    for ll=1:nummom+1
        l=ll-1;
        ld2=floor(l/2);
        evenl=double(mod(l,2)==0);
        % Calculate numerical coefficients Am and Bi in Dave, double-sums for moments
        if l==0
            idel=1;
            am(1:ntrm+1)=2*recip(2*(0:1:ntrm)+1);
            bi(1)=1;
        elseif (evenl==1) && (l~=0)
            idel=1;
            for m=ld2:ntrm
                mm=m+1;
                am(mm)=(1+recip(2*m-l+1))*am(mm);
            end
            for i=0:ld2-1
                ii=i+1;
                bi(ii)=(1-recip(l-2*i))*bi(ii);
            end
            bi(ld2+1)=(2-recip(l))*bi(ld2);
        else
            idel=2;
            for m=ld2:ntrm
                mm=m+1;
                am(mm)=(1-recip(2*m+l+2))*am(mm);
            end
            for i=0:1:ld2
                ii=i+1;
                bi(ii)=(1-recip(l+2*i+1))*bi(ii);
            end
        end
        % establish upper limits for sums and incorporate factor DEL into Bi
        mmax=ntrm-idel;
        if ipolzn>=0, mmax=mmax+1; end
        imax=min([ld2 mmax-ld2]);
        if imax<0
            return;
        end
        bidel(1:imax+1)=bi(1:imax+1);
        if evenl==1, bidel(1)=0.5*bidel(1); end
        if ipolzn==0 % unpolarized phase function: (M1+M2)/2
            for i=0:1:imax
                ii=i+1;
                tsum=0;
                for m=ld2:mmax-i
                    mm=m+1;
                    tsum=tsum+am(mm)*(real(c(m-i+1)*conj(c(m+i+idel)))+real(d(m-i+1)*conj(d(m+i+idel))));
                end
                pmom(ll,1)=pmom(ll,1)+bidel(ii)*tsum;
            end
            pmom(ll,1)=0.5*pmom(ll,1);
        else
            if calcmo(1)==1 % M1 (Sekera R1 if ipolzn<0): sums of C*conj(C)
                for ii=1:imax+1
                    i=ii-1;
                    tsum=0;
                    for mm=ld2+1:mmax-i+1
                        m=mm-1;
                        tsum=tsum+am(mm)*real(c(m-i+1)*conj(c(m+i+idel)));
                    end
                    pmom(ll,1)=pmom(ll,1)+bidel(ii)*tsum;
                end
            end
            if calcmo(2)==1 % M2 (R2): sums of D*conj(D)
                for ii=1:imax+1
                    i=ii-1;
                    tsum=0;
                    for mm=ld2+1:mmax-i+1
                        m=mm-1;
                        tsum=tsum+am(mm)*real(d(m-i+1)*conj(d(m+i+idel)));
                    end
                    pmom(ll,2)=pmom(ll,2)+bidel(ii)*tsum;
                end
            end
            if calcmo(3)==1 % S21 (R3): real part of C*conj(D) sums
                for ii=1:imax+1
                    i=ii-1;
                    tsum=0;
                    for mm=ld2+1:mmax-i+1
                        m=mm-1;
                        tsum=tsum+am(mm)*(real(c(m-i+1)*conj(d(m+i+idel)))+real(c(m+i+idel)*conj(d(m-i+1))));
                    end
                    pmom(ll,3)=pmom(ll,3)+bidel(ii)*tsum;
                end
                pmom(ll,3)=0.5*pmom(ll,3);
            end
            if calcmo(4)==1 % D21 (R4): imaginary part of C*conj(D) sums
                for ii=1:imax+1
                    i=ii-1;
                    tsum=0;
                    for mm=ld2+1:mmax-i+1
                        m=mm-1;
                        tsum=tsum+am(mm)*(imag(c(m-i+1)*conj(d(m+i+idel)))+imag(c(m+i+idel)*conj(d(m-i+1))));
                    end
                    pmom(ll,4)=pmom(ll,4)+bidel(ii)*tsum;
                end
                pmom(ll,4)=-0.5*pmom(ll,4);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmom=lpco1t(nmom,ipolzn,~,calcmo,a,b)

% function pmom=lpco1t(nmom,ipolzn,momdim,calcmo,a,b)
% Calculate Legendre polynomial expansion coefficients (also called moments)
%  for phase quantities in special case where no. terms in Mie series = 1
% Subroutine used in MIEV0 (Fortran routine LPCO1T), called from lpcoef
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 1 May 2012
% Bug fixes 8 October 2026 (with Claude): pmom is now preallocated as
%  (nmom+1) x (highest requested phase quantity), so its shape no longer
%  depends on momdim; Mueller loop limited to moments 0-2 as in the
%  Fortran (higher moments are zero for a 1-term series).
%
% Inputs:  nmom, ipolzn, calcmo - as in miev0; momdim - unused
%          a, b - Mie coefficients a_n, b_n (only a(1), b(1) used)
% Output:  pmom - (nmom+1) x k Legendre moments, row 1 = moment 0
% Functions used:  sq

% preallocation
pmom=zeros(nmom+1,max([1 find(calcmo,1,'last')]));

a1sq=sq(a(1));
b1sq=sq(b(1));
a1b1c=a(1)*conj(b(1));

if ipolzn<0
    if calcmo(1)==1, pmom(1,1)=2.25*b1sq; end
    if calcmo(2)==1, pmom(1,2)=2.25*a1sq; end
    if calcmo(3)==1, pmom(1,3)=2.25*real(a1b1c); end
    if calcmo(4)==1, pmom(1,4)=2.25*imag(a1b1c); end
else
    nummom=min([nmom 2]);
    for ll=1:nummom+1
        l=ll-1;
        if ipolzn==0
            if l==0, pmom(ll,1)=1.5*(a1sq+b1sq); end
            if l==1, pmom(ll,1)=1.5*real(a1b1c); end
            if l==2, pmom(ll,1)=0.15*(a1sq+b1sq); end
            continue
        end
        if calcmo(1)==1
            if l==0, pmom(ll,1)=2.25*(a1sq+b1sq/3); end
            if l==1, pmom(ll,1)=1.5*real(a1b1c); end
            if l==2, pmom(ll,1)=0.3*b1sq; end
        end
        if calcmo(2)==1
            if l==0, pmom(ll,2)=2.25*(b1sq+a1sq/3); end
            if l==1, pmom(ll,2)=1.5*real(a1b1c); end
            if l==2, pmom(ll,2)=0.3*a1sq; end
        end
        if calcmo(3)==1
            if l==0, pmom(ll,3)=3*real(a1b1c); end
            if l==1, pmom(ll,3)=0.75*(a1sq+b1sq); end
            if l==2, pmom(ll,3)=0.3*real(a1b1c); end
        end
        if calcmo(4)==1
            if l==0, pmom(ll,4)=-1.5*imag(a1b1c); end
            if l==1, pmom(ll,4)=0; end
            if l==2, pmom(ll,4)=0.3*imag(a1b1c); end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmom=lpco2t(nmom,ipolzn,~,calcmo,a,b)

% function pmom=lpco2t(nmom,ipolzn,momdim,calcmo,a,b)
% Calculate Legendre polynomial expansion coefficients (also called moments)
%  for phase quantities in special case where no. terms in Mie series = 2
% Subroutine used in MIEV0 (Fortran routine LPCO2T), called from lpcoef
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 1 May 2012
% Bug fixes 8 October 2026 (with Claude): pmom is now preallocated as
%  (nmom+1) x (highest requested phase quantity), so its shape no longer
%  depends on momdim.
%
% Inputs:  nmom, ipolzn, calcmo - as in miev0; momdim - unused
%          a, b - Mie coefficients a_n, b_n (a(1:2), b(1:2) used)
% Output:  pmom - (nmom+1) x k Legendre moments, row 1 = moment 0
%          (at most moments 0-4 are nonzero for a 2-term series)
% Functions used:  sq

% preallocation
pmom=zeros(nmom+1,max([1 find(calcmo,1,'last')]));

ca=3*a(1)-5*b(2);
cat=3*b(1)-5*a(2);
cac=conj(ca);
a2sq=sq(a(2));
b2sq=sq(b(2));
a2c=conj(a(2));
b2c=conj(b(2));

if ipolzn<0
    nummom=min([nmom 2]);
    for ll=1:nummom+1 % loop over Sekera moments
        l=ll-1;
        if calcmo(1)==1
            if l==0, pmom(ll,1)=0.25*(sq(cat)+(100/3)*b2sq); end
            if l==1, pmom(ll,1)=(5/3)*real(cat*b2c); end
            if l==2, pmom(ll,1)=(10/3)*b2sq; end
        end
        if calcmo(2)==1
            if l==0, pmom(ll,2)=0.25*(sq(ca)+(100/3)*a2sq); end
            if l==1, pmom(ll,2)=(5/3)*real(ca*a2c); end
            if l==2, pmom(ll,2)=(10/3)*a2sq; end
        end
        if calcmo(3)==1
            if l==0, pmom(ll,3)=0.25*real(cat*cac+(100/3)*b(2)*a2c); end
            if l==1, pmom(ll,3)=(5/6)*real(b(2)*cac+cat*a2c); end
            if l==2, pmom(ll,3)=(10/3)*real(b(2)*a2c); end
        end
        if calcmo(4)==1
            if l==0, pmom(ll,4)=-0.25*imag(cat*cac+(100/3)*b(2)*a2c); end
            if l==1, pmom(ll,4)=(-5/6)*imag(b(2)*cac+cat*a2c); end
            if l==2, pmom(ll,4)=(-10/3)*imag(b(2)*a2c); end
        end
    end
else
    cb=3*b(1)+5*a(2);
    cbt=3*a(1)+5*b(2);
    cbc=conj(cb);
    cg=(cbc*cbt+10*(cac*a(2)+b2c*cat))/3;
    ch=2*(cbc*a(2)+b2c*cbt);
    nummom=min([nmom 4]);
    for ll=1:nummom+1 % loop over Mueller moments
        l=ll-1;
        if (ipolzn==0) || (calcmo(1)==1)
            if l==0, pm1=0.25*sq(ca)+sq(cb)/12+(5/3)*real(ca*b2c)+5*b2sq; end
            if l==1, pm1=real(cb*(cac/6+b2c)); end
            if l==2, pm1=sq(cb)/30+(20/7)*b2sq+(2/3)*real(ca*b2c); end
            if l==3, pm1=(2/7)*real(cb*b2c); end
            if l==4, pm1=(40/63)*b2sq; end
            if calcmo(1)==1, pmom(ll,1)=pm1; end
        end
        if (ipolzn==0) || (calcmo(2)==1)
            if l==0, pm2=0.25*sq(cat)+sq(cbt)/12+(5/3)*real(cat*a2c)+5*a2sq; end
            if l==1, pm2=real(cbt*(conj(cat)/6+a2c)); end
            if l==2, pm2=sq(cbt)/30+(20/7)*a2sq+(2/3)*real(cat*a2c); end
            if l==3, pm2=(2/7)*real(cbt*a2c); end
            if l==4, pm2=(40/63)*a2sq; end
            if calcmo(2)==1, pmom(ll,2)=pm2; end
        end
        if ipolzn==0
            pmom(ll,1)=0.5*(pm1+pm2);
            continue
        end
        if calcmo(3)==1
            if l==0, pmom(ll,3)=0.25*real(cac*cat+cg+20*b2c*a(2)); end
            if l==1, pmom(ll,3)=real(cac*cbt+cbc*cat+3*ch)/12; end
            if l==2, pmom(ll,3)=0.1*real(cg+(200/7)*b2c*a(2)); end
            if l==3, pmom(ll,3)=real(ch)/14; end
            if l==4, pmom(ll,3)=(40/63)*real(b2c*a(2)); end
        end
        if calcmo(4)==1
            if l==0, pmom(ll,4)=0.25*imag(cac*cat+cg+20*b2c*a(2)); end
            if l==1, pmom(ll,4)=imag(cac*cbt+cbc*cat+3*ch)/12; end
            if l==2, pmom(ll,4)=0.1*imag(cg+(200/7)*b2c*a(2)); end
            if l==3, pmom(ll,4)=imag(ch)/14; end
            if l==4, pmom(ll,4)=(40/63)*imag(b2c*a(2)); end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output=sq(ctmp)

% function output=sq(ctmp)
% Squared modulus of a complex number, |ctmp|^2 (elementwise)
% Subroutine used in MIEV0 (statement function SQ in the Fortran)
%
% Original FORTRAN by Warren Wiscombe
% Matlabized by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% 1 May 2012
%
% Inputs:  ctmp
% Output:  output

output=real(ctmp).^2+imag(ctmp).^2;
