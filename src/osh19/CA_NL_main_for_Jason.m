% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Main script to run to solve the nonlinear FARE model
% with rotation.  See report for important information regarding how
% variables are defined and methods used

% clear all

disp('------------------------------------------------');
disp('------------------------------------------------');
disp('Running CA_NLmodel/CA_NL_main_for_Jason');
disp('------------------------------------------------');
disp('------------------------------------------------');

% Make sure script is run from the right folder
% cd /Users/hrogrosky/Documents/CA_NLmodel;

% Run script to initialize all values and parameters, including IC
CA_NL_init;
disp('------------------------------------------------');
disp('Marching forward in time');
disp('------------------------------------------------');

% Setup saved solution matrices
% This version keeps track of baroclinic winds at all levels (except the
% top level) and barotropic winds, rather than total winds at each level.
usol=zeros(numy,numx,numlevs+1,tnumtosave);
vsol=zeros(numy,numx,numlevs+1,tnumtosave);
zsol=zeros(numy,numx,1,tnumtosave);
psol=zeros(numy,numx,numlevs+1,tnumtosave);
wsol=zeros(numy,numx,numlevs+1,tnumtosave);
thsol=zeros(numy,numx,numlevs+1,tnumtosave);
qsol=zeros(numy,numx,numlevs+1,tnumtosave);
qsolanom=zeros(numy,numx,numlevs+1,tnumtosave);
thsolanom=zeros(numy,numx,numlevs+1,tnumtosave);

% But we'll also save *total* winds and pressure as we go along
utotsol=zeros(numy,numx,numlevs+1,tnumtosave);
vtotsol=zeros(numy,numx,numlevs+1,tnumtosave);

% Initialize a few parameters - GENERALLY SHOULDN'T CHANGE THESE!
runsolver=1;
tcount=0;
timestepcount=0; numtimestepstocalcavgtime=100;
solnum=0;
dthalver=0;
tnumtosave_v2=tnumtosave;
timepertimestep=zeros(1,numtimestepstocalcavgtime);

% Solve model equations
if runsolver==1
    % Step forward in time
    for tn=tt %(1:10)
        timestepcount=timestepcount+1;
        if timestepcount<numtimestepstocalcavgtime
            tic
        end
        if timestepcount==numtimestepstocalcavgtime
            disp('------------------------------------------------');
            disp(strcat(['Average computing time per simulation day=', num2str(sum(timepertimestep)/numtimestepstocalcavgtime*saveevery),' s']));
            disp(strcat(['Expected total computing time=', num2str(sum(timepertimestep)/numtimestepstocalcavgtime*saveevery*(tnumtosave-1)/60),' min']));
            disp('------------------------------------------------');
        end
        % Periodically save solutions
        if mod(tcount,saveevery)==0
            solnum=solnum+1; %round(tn/tsdt+1);
            usol(:,:,:,solnum)=u;
            vsol(:,:,:,solnum)=v;
            zsol(:,:,1,solnum)=z;
            wsol(:,:,:,solnum)=w;
            psol(:,:,:,solnum)=p;
            thsol(:,:,:,solnum)=th;
            thsolanom(:,:,:,solnum)=thanom; %-thetaetildemat;
            qsol(:,:,:,solnum)=q;
            qsolanom(:,:,:,solnum)=qanom;
            [utot,vtot]=CA_NL_calc_total_winds( u,v,z,beta,YY,scalex,scaley,dafrac );
            utotsol(:,:,:,solnum)=utot;
            vtotsol(:,:,:,solnum)=vtot;
            disp(strcat(['Day ',int2str(solnum-1),'/',int2str(tnumtosave_v2-1),', Max u is ',num2str(1000*max(max(max(abs(utot))))),' m/s']));
            % Check that you have real-valued data
            if isnan(u(1,1,numlevs))
                disp('Warning:  NaNs!!');
            end
        end
        
        
        % March forward in time
        theta=th; Y=YY;
        CA_NL_2BCCA_RK4_method;
        u=utemp; v=vtemp; z=ztemp; th=thetatemp; q=qtemp; w=wtemp; p=ptemp; qanom=qanomtemp; thanom=thanomtemp;
        
%         % Print output of the first step to base workspace.
%         assignin('base', 'uOrig', u);
%         assignin('base', 'vOrig', v);
%         assignin('base', 'zOrig', z);
%         assignin('base', 'tOrig', thanom);
%         assignin('base', 'qOrig', qanom);
%         disp('ORIG: Output first step to base workspace!');
%         pause(inf);
        
        tcount=tcount+2^(-dthalver);
        if timestepcount<numtimestepstocalcavgtime
            timepertimestep(timestepcount)=toc;
        end
    end
end

% Save the results
disp('Saving data');
CA_NL_save;
disp('------------------------------------------------');
disp('Finished');
disp('------------------------------------------------');

% % Plotting scripts
% CA_NL_make_horcs_figs;
% CA_NL_make_vertlong_figs;
% CA_NL_make_timelong_figs;
% CA_NL_make growthrate_figs;
