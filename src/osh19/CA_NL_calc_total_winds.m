% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Script that takes in *baroclinic* winds and 
% *barotropic* relative vorticity, and returns *total* winds at levels 
% j=1,...,J

function [ utot,vtot ] = CA_NL_calc_total_winds( uBC,vBC,zeta,beta,Y,scalex,scaley,dafrac )

    [numy,numx,numz]=size(uBC);
    numlevs=numz-1;
    
    % Calculate BT variables
    [ubar,vbar,~,~,~,~]=CA_NL_calc_pbar(uBC,vBC,zeta,beta,Y,scalex,scaley,dafrac);
    
    % Calculate *total* winds at each level j=1,...,J-1
    utot=zeros(numy,numx,numlevs+1);
    vtot=zeros(numy,numx,numlevs+1);
    for icount=2:numlevs
        utot(:,:,icount)=uBC(:,:,icount)+ubar;
        vtot(:,:,icount)=vBC(:,:,icount)+vbar;
    end 

    % Calculate *total* winds at top level, j=J
    utot(:,:,numlevs+1)=numlevs*ubar-squeeze(sum(utot(:,:,2:numlevs),3));
    vtot(:,:,numlevs+1)=numlevs*vbar-squeeze(sum(vtot(:,:,2:numlevs),3));

end

