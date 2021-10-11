% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  RK4 method for time-marching when solving
% nonlinear FARE model.  See report for additional information.

% Predictor step 1 (Y2)
utemp=u; vtemp=v; ztemp=z; qtemp=q; ptemp=p; wtemp=w; qanomtemp=qanom; thanomtemp=thanom; % thtemp=th??

% CA_NL_RHS_2BCCAsw;
CA_NL_RHS_2BCCA;
% CA_NL_RHS_2BCCA_sw_12_18;
uRHS1=uRHS; vRHS1=vRHS; zRHS1=zRHS; thetaRHS1=thetaRHS; qRHS1=qRHS;

utemp=u+dt/2*uRHS1;      
vtemp=v+dt/2*vRHS1;
ztemp=z+dt/2*zRHS1;
% thetatemp=th+dt/2*thetaRHS1;
% qtemp=q+dt/2*qRHS1;
qanomtemp=qanom+dt/2*qRHS1;
thanomtemp=thanom+dt/2*thetaRHS1;
[wtemp,ptemp]=CA_NL_calc_wp(utemp,vtemp,g,theta0,dz,scalex,scaley,qbgmat,dafrac,qanomtemp,thanomtemp);

% Save w, p params to workspace.
% assignin('base', 'uclinOrig', utemp);
% assignin('base', 'vclinOrig', vtemp);
% assignin('base', 'thanomOrig', thanomtemp);
% assignin('base', 'dzOrig', dz);
% assignin('base', 'gOrig', g);
% assignin('base', 'theta0Orig', theta0);
% assignin('base', 'scalexOrig', scalex);
% assignin('base', 'scaleyOrig', scaley);
% assignin('base', 'wOrig', wtemp);
% assignin('base', 'pOrig', ptemp);
% disp('Orig: Output w, p params to base workspace!');
% pause(inf);

% Predictor step 2 (Y3)
% CA_NL_RHS_2BCCAsw;
CA_NL_RHS_2BCCA;
% CA_NL_RHS_2BCCA_sw_12_18;
uRHS2=uRHS; vRHS2=vRHS; zRHS2=zRHS; thetaRHS2=thetaRHS; qRHS2=qRHS;



utemp=u+dt/2*uRHS2;      
vtemp=v+dt/2*vRHS2;
ztemp=z+dt/2*zRHS2;
% thetatemp=th+dt/2*thetaRHS2;
% qtemp=q+dt/2*qRHS2;
qanomtemp=qanom+dt/2*qRHS2;
thanomtemp=thanom+dt/2*thetaRHS2;
[wtemp,ptemp]=CA_NL_calc_wp(utemp,vtemp,g,theta0,dz,scalex,scaley,qbgmat,dafrac,qanomtemp,thanomtemp);


% Predictor step 3 (Y4)
% CA_NL_RHS_2BCCAsw;
CA_NL_RHS_2BCCA;
% CA_NL_RHS_2BCCA_sw_12_18;
uRHS3=uRHS; vRHS3=vRHS; zRHS3=zRHS; thetaRHS3=thetaRHS; qRHS3=qRHS;

utemp=u+dt*uRHS3;      
vtemp=v+dt*vRHS3;
ztemp=z+dt*zRHS3;
% thetatemp=th+dt*thetaRHS3;
% qtemp=q+dt*qRHS3;
qanomtemp=qanom+dt*qRHS3;
thanomtemp=thanom+dt*thetaRHS3;

% % Save RHS to base workspace.
% assignin('base', 'uOrig', uRHS3);
% assignin('base', 'vOrig', vRHS3);
% assignin('base', 'zOrig', zRHS3);
% assignin('base', 'tOrig', thetaRHS3);
% assignin('base', 'qOrig', qRHS3);
% disp('ORIG: Output RK4 RHS to base workspace!');
% pause(inf);

[wtemp,ptemp]=CA_NL_calc_wp(utemp,vtemp,g,theta0,dz,scalex,scaley,qbgmat,dafrac,qanomtemp,thanomtemp);

% % Save RHS to base workspace.
% assignin('base', 'wOrig', wtemp);
% assignin('base', 'pOrig', ptemp);
% disp('ORIG: Output w, p to base workspace!');
% pause(inf);

% Predictor step 4 (f(Y4))
% CA_NL_RHS_2BCCAsw;
CA_NL_RHS_2BCCA;
% CA_NL_RHS_2BCCA_sw_12_18;

% % Save RHS to base workspace.
% assignin('base', 'uOrig', uRHS);
% assignin('base', 'vOrig', vRHS);
% assignin('base', 'zOrig', zRHS);
% assignin('base', 'tOrig', thetaRHS);
% assignin('base', 'qOrig', qRHS);
% disp('ORIG: Output RK4 RHS to base workspace!');
% pause(inf);

% March forward in time
utemp=u+dt/6*(uRHS1+2*uRHS2+2*uRHS3+uRHS);
vtemp=v+dt/6*(vRHS1+2*vRHS2+2*vRHS3+vRHS);
ztemp=z+dt/6*(zRHS1+2*zRHS2+2*zRHS3+zRHS);
% thetatemp=th+dt/6*(thetaRHS1+2*thetaRHS2+2*thetaRHS3+thetaRHS);
% qtemp=q+dt/6*(qRHS1+2*qRHS2+2*qRHS3+qRHS);
qanomtemp=qanom+dt/6*(qRHS1+2*qRHS2+2*qRHS3+qRHS);
qtemp=qbgmat+qanomtemp;
thanomtemp=thanom+dt/6*(thetaRHS1+2*thetaRHS2+2*thetaRHS3+thetaRHS);
thetatemp=thetatildemat+thanomtemp;
[wtemp,ptemp]=CA_NL_calc_wp(utemp,vtemp,g,theta0,dz,scalex,scaley,qbgmat,dafrac,qanomtemp,thanomtemp);

clear uRHS1 uRHS2 uRHS3 uRHS
clear vRHS1 vRHS2 vRHS3 vRHS
clear zRHS1 zRHS2 zRHS3 zRHS
% clear thetaRHS1 thetaRHS2 thetaRHS3 thetaRHS
clear qRHS1 qRHS2 qRHS3 qRHS
