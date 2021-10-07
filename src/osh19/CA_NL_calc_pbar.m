% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Calculate barotropic pressure from barotropic winds 
% by solving Poisson problem (spectrally)

function [ubar,vbar,u2bar,uvbar,v2bar,pbar]=CA_NL_calc_pbar(u,v,z,beta,Y,scalex,scaley,dafrac)
    [numy,numx,numz]=size(u);
    numlevs=numz-1;

    wavenumbervecy=transpose([0:numy/2 -numy/2+1:-1])*scaley;
    wavenumbervecx=[0:numx/2 -numx/2+1:-1]*scalex;
    [kk,ll]=meshgrid(wavenumbervecx,wavenumbervecy);
    wavenumbermat=-(ll.^2+kk.^2);
    wavenumbermat(1,1)=10^10;  % THIS NUMBER SHOULD BE IRRELEVANT, AS LONG AS IT'S NONZERO, RIGHT?  
    % (THE RHS OF THE PBAR EQUATION SHOULD SATISFY A ZERO MEAN CONDITION SO THAT IT SHOULDN'T MATTER) 

    % Calculate *barotropic* winds from barotropic relative vorticity
    psibar=ifft2(fft2(z)./wavenumbermat);
    psibar=dealias(psibar,dafrac);
    ubar=-D1(psibar,'y',scaley);
    vbar=D1(psibar,'x',scalex);

    % Calculate *total* winds at each level j=1,...,J-1
    utot=zeros(size(u));
    vtot=zeros(size(u));
    for icount=2:numlevs
        utot(:,:,icount)=u(:,:,icount)+ubar;
        vtot(:,:,icount)=v(:,:,icount)+vbar;
    end 

    % Calculate *top-level total* winds from barotropic winds
    utoptot=numlevs*ubar-squeeze(sum(utot(:,:,2:end),3));
    vtoptot=numlevs*vbar-squeeze(sum(vtot(:,:,2:end),3));

    % Calculate u2bar, uvbar, v2bar
    u2bar=1/numlevs*(sum(utot(:,:,2:end).*utot(:,:,2:end),3)+utoptot.*utoptot);
    uvbar=1/numlevs*(sum(utot(:,:,2:end).*vtot(:,:,2:end),3)+utoptot.*vtoptot);
    v2bar=1/numlevs*(sum(vtot(:,:,2:end).*vtot(:,:,2:end),3)+vtoptot.*vtoptot);

    % Calculate *barotropic* pressure - calculating derivative of beta*y
    % numerically
    pbarRHS=-D2(u2bar,'x','x',scalex,scaley)-2*D2(uvbar,'x','y',scalex,scaley)-D2(v2bar,'y','y',scalex,scaley)...
        +beta*Y.*D1(vbar,'x',scalex)-beta*D1(Y.*ubar,'y',scaley);

    pbar=ifft2(fft2(pbarRHS)./wavenumbermat);
    pbar=dealias(pbar,dafrac);
end
