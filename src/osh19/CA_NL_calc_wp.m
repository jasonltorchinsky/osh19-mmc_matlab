% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Calculate w and p from u, v, theta, qv for nonlinear
% FARE model

function [w,p]=CA_NL_calc_wp(u,v,g,theta0,dz,scalex,scaley,qbgtildemat,dafrac,qanom,thanom)
    % NOTE:  Input variables are all *baroclinic*, except for zeta 
    % (barotropic relative vorticity)

    [~,~,numz]=size(u);
    numlevs=numz-1;
    p=zeros(size(u));
    w=zeros(size(u));
    
%     [~,~,~,~,~,pbar]=CA_NL_calc_pbar(u,v,z,beta,Y,scalex,scaley,dafrac);
%     p(:,:,2)=pbar;
    
    % Calculate pressure
    for kcount=1:numlevs-1
        p(:,:,2)=p(:,:,2)-g*dz/numlevs*(numlevs-kcount)*(1/theta0*(thanom(:,:,kcount+1))); %...
%             -(qanom(:,:,kcount+1)-qbgtildemat(:,:,kcount+1))); %+(eps0-Lv/cp/theta0)*(qv(:,:,kcount+1)-qvsmat(:,:,kcount+1))
    end
    for kcount=3:numlevs+1
        p(:,:,kcount)=p(:,:,kcount-1)+g*dz*(1/theta0*(thanom(:,:,kcount-1))); %...
%             -(qanom(:,:,kcount-1)-qbgtildemat(:,:,kcount-1))); %+(eps0-Lv/cp/theta0)*(qv(:,:,kcount-1)-qvsmat(:,:,kcount-1))
    end

    % Calculate vertical winds
    for jcount=2:numlevs
        w(:,:,jcount)=w(:,:,jcount-1)-dz*(D1(u(:,:,jcount),'x',scalex)+D1(v(:,:,jcount),'y',scaley));
    end
end

