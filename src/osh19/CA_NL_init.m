% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Define all parameter values and initialize all
% variables for solving nonlinear FARE model

% close all
figfontsize=17;
printyes = 1;

%%%%%%%% CHANGE THE PLACE AND NAME THE INITIAL CONDITIONS ARE SAVED!! %%%%%

%%%%%%%% USER-SPECIFIED PARAMETERS %%%%%%%%

%%%%%%%%%% CURRENTLY IN CASE 5 %%%%%%%%%%%

numx=32;                  % Number of x gridpoints
numy=32;                  % Number of y gridpoints
numlevs=8;                % Number of nontrivial levels for (u,v,p)  
tend=100;                 % Number of simulation days
pY=6000;                  % Domain width in y (km)                                       DEFAULT = 6000
Bvsdim=-0.00134;%-0.00134;          % (kg kg^{-1} km^{-1}) From HSS15 - range of -1 to -1.5
tau_u=25;%25;                 % Damping timescale in days - u, v, zeta
tau_theta=25;%25;             % Damping timescale in days - q
tautop=1;                 % days    
taubottom=2/24;%2/24;           % days    
btop=0.8*76;%0.8*76;              % bmid=0.8 in SH17 - multiply by ~3 to make dimensional??  
bbottom=0.1*76;%0.1*76;           % blow=0.1 in SH17 - multiply by ~3 to make dimensional??
Dvtop=0.0001;%0.0001;       
Dvbottom=0.0001;%0.0001;
QBGvecversion=2;
tauvecversion=0;
bvecversion=0;
Dvversion=0;
use_2BCConvAdj=2;
fnameIC='Stan';
expscale=120000;%120000; 
Lscale=2000;%2000;
BvsYvecamp=0.25;%0.25; 

IC_type=3;                % Choose your initial condition
                          % 1 = Single mode, single wavenumber
                          % 2 = Single mode, multiple wavenumbers
                          % 3 = Multiple modes, multiple wavenumbers
                          % 4 = Load results from previous nonlinear simulation (NOT SETUP CORRECTLY YET!)
                          % 5 = Load previously computed initial condition

%%%%%%%% END USER-SPECIFIED PARAMETERS %%%%%%%%
%%%%%%%% BUT SEE BELOW FOR INITIAL CONDITIONS CHOICE %%%%%%%%


%%%%%%%% 'FIXED' PARAMETERS %%%%%%%%
H=16;                     % Height of troposphere (km)
dz=H/numlevs;             % km
beta=2.3*10^(-8);         % 1/(km*s)
g=9.8*10^(-3);            % km/s^2
Bdim=3;                   % K km^{-1} From HSS15 - they use a range of 2-4
cp = 1000;                % J kg^{-1} K^-1
theta0=300;               % K
Lv = 2.5*10^6;            % J kg^{-1}
eps0=0.6;                 
c=sqrt(g*Bdim/theta0)*H/pi; % Wave speed of the system

QBCtop=0;
QBCbottom=-Bvsdim*H;
zvec=transpose(dz:dz:H-dz);
sinzvec=transpose(sin(zvec*pi/H));
rsinzvec=reshape(sinzvec,[1,1,length(sinzvec)]);
sinzmat=repmat(rsinzvec,[numy,numx,1]);
sin2zvec=transpose(sin(2*zvec*pi/H));
rsin2zvec=reshape(sin2zvec,[1,1,length(sin2zvec)]);
sin2zmat=repmat(rsin2zvec,[numy,numx,1]);
atemp=QBCbottom/(1-exp(-H/expscale));
btemp=-atemp*exp(-H/expscale);
QBGvec=atemp*exp(-zvec/expscale)+btemp;
QBGvecwithBC=[QBCbottom; QBGvec; QBCtop];
Bvsvec=1/(2*dz)*(QBGvecwithBC(3:end)-QBGvecwithBC(1:end-2));


tau_upper=tautop;
tau_lower=taubottom;
if tauvecversion==0         % Linear ThetaBGvec (old version);
    tauvec=tautop+(tautop-taubottom)*(zvec-H)/H;
end
if bvecversion==0           % Linear bvec
    bvec=btop+(btop-bbottom)*(zvec-H)/H;
end
if Dvversion==0             % Linear bvec
    Dvvec=Dvtop+(Dvtop-Dvbottom)*(zvec-H)/H;
end

% Constants (length/temporal scales)
Pe = 40000;               % km - circumference of the earth. 
scalex=2*pi/Pe;           
scaley=2*pi/(2*pY);

% Convert damping timescales from days to seconds
tauvec=tauvec*24*3600;      % seconds
tau_upper=tau_upper*24*3600;
tau_lower=tau_lower*24*3600;
tau_u =tau_u *24*3600;      % seconds
tau_theta=tau_theta*24*3600; % seconds

% Some other parameters and some rescaling of parameters to make units
% consistent
alphabar=H/pi*Bdim;       % K
Q=cp*alphabar/Lv;
Fscale=H/pi*(Bdim+Lv/cp*Bvsdim-theta0*Bvsdim);
Ftilde=Fscale/alphabar;
G=H/(pi*Q)*Bvsdim/Ftilde;
Fsat=Ftilde;
cmoist=c*sqrt(Ftilde);

% Constants (length/temporal scales)
L = sqrt(c/beta);         % km

% Dealiasing parameter
dafrac=2/3;
%%%%%%%% END 'FIXED' PARAMETERS %%%%%%%%




%%%%%%%% CHOOSE INITIAL CONDITION TYPE %%%%%%%%
% Some amplitude and phase shift multipliers
IC_amps_by_k=[0.7547 0.2760 0.6797 0.6551 0.1626 0.1190 0.4984 0.9597 0.3404 0.5853 0.2238 0.7513 0.2551 0.5060 0.6991 0.8909 0.9593 0.5472 0.1386 0.1493];
IC_ps_by_k=[1.0180 5.2824 1.5977 5.1163 1.5301 5.8387 2.1990 1.2352 1.5776 3.8707 2.9738 2.2095 5.2203 3.6773 3.4540 5.7629 1.7960 4.7576 4.7358 2.3904];

if IC_type==1              % Use for single wavenumber, single mode
    % Pick an initial condition from the linear solutions
    % Number of the mode you want to use - (sorted by growthrates, largest to smallest)
    num_damped = 1;        % Which mode do you want?  1=most unstable, 6=most damped.  DEFAULT = 1!!  
    wave_number = 2;       % What zonal wavenumber do you want?  DEFAULT = 2
    amp_factor=10^(-1);     % How large do you want the initial condition?  
    
elseif IC_type==2          % Use for combining multiple wavenumbers, single mode
    % Pick an initial condition from the linear solutions
    % Number of the mode you want to use - (sorted by growthrates, largest to smallest)
    num_damped = 1;        % Which mode do you want?  1=most unstable, 6=most damped.  DEFAULT = 1!!  
    st_wave_number = 1;    % What zonal wavenumber do you want?  DEFAULT = 2
    amp_factor=3*10^(-1);    % How large do you want the initial condition?  
    numwaves_in_IC=3;      % How many wavenumbers do you want to add together in the initial condition?
    
elseif IC_type==3          % Use for combining multiple modes and multiple wavenumbers
    % Pick an initial condition from the linear solutions
    % Number of the mode you want to use - (sorted by growthrates, largest to smallest)
    num_damped_vec = [1 2 3 4]; % 2 3 4];        % Which modes do you want?  1=most unstable, 6=most damped.  DEFAULT = 1!!  
    st_wave_number = 1;    % What zonal wavenumber do you want?  DEFAULT = 2
    amp_factor=1*10^(0);    % How large do you want the initial condition?  
    numwaves_in_IC=3;      % How many wavenumbers do you want to add together in the initial condition?
    disp('------------------------------------------------');
    disp(strcat(['Creating IC using ',int2str(length(num_damped_vec)),' modes with k=',int2str(st_wave_number),'-',int2str(st_wave_number+numwaves_in_IC-1)]));
    disp('------------------------------------------------');
    
elseif IC_type==4          % Use for restarting from previous nonlinear run
%     FARE_load_previous_run;

elseif IC_type==5          % Load a previously-calculated initial condition
    CA_NL_load_data;

end
%%%%%%%% END CHOOSE INITIAL CONDITION TYPE %%%%%%%%

%%%%%%%% SET UP SPACE AND TIME GRIDS (AND REMAINING PARAMETERS) %%%%%%%%
% Setup grid to solve on
dx=Pe/numx;
xx=0:dx:Pe-dx;
dy=2*pY/numy;
yy=-pY+dy/2:dy:pY-dy/2;
zzU=-dz/2:dz:H-dz/2;
zzW=0:dz:H;
[XX YY]=meshgrid(xx,yy);
[XX3U YY3U ZZ3U]=meshgrid(xx,yy,zzU);
[XX3W YY3W ZZ3W]=meshgrid(xx,yy,zzW);

% Create amplitude function for Bvs (possibly with y-dependence)
BvsYvec=ones(size(yy))-BvsYvecamp*(1-exp(-(yy/Lscale).^2/2));

% Setup time grid
tend=tend*3600*24;

% Choose time interval for marching
CFL_constant=0.15/40; % 0.15?  Not necessary?  
dt_x=CFL_constant*dx/cmoist; % Calculates dt in seconds
dt_temp=dt_x;

% Round down to nearest fraction of a day
dt_day=dt_temp/3600/24;
dt_day_round=ceil(1/dt_day);
dt=1/dt_day_round*3600*24;
dtstart=dt;

t0=0;
tt=0:dt:tend;
tnum=length(tt)-1;

% Choose time interval for saving solution
saveeverynthday=1;
saveevery=round(saveeverynthday/(dt/3600/24));
tsdt=dt*saveevery;
tts=tt(1:saveevery:end);
tnumtosave=length(tts);

% Set up moisture flux term - set to zeros (or other constant) to test
% linear theory results
F=zeros(size(zzW));

% Set up background vertical moisture gradient - set to linear profile to
% test linear theory results
qbg=QBGvecwithBC;
[~,~,qbgmat]=meshgrid(xx,yy,qbg);
for xcount=1:numx
    for zcount=1:numlevs+1
        qbgmat(:,xcount,zcount)=squeeze(qbgmat(:,xcount,zcount)).*transpose(BvsYvec);  % CHANGE 7/18/17!!
    end
end
qbg=qbgmat;

% Set up background potential temperature vertical gradient - set to linear
% profile to test linear theory results
thetatilde=Bdim*(zzW);
[~,~,thetatildemat]=meshgrid(xx,yy,thetatilde); % CHANGE 7/18/17!!

u=zeros(numy,numx,numlevs+1);
v=zeros(numy,numx,numlevs+1);
z=zeros(numy,numx,1);
thanom=zeros(numy,numx,numlevs+1);
qanom=zeros(numy,numx,numlevs+1);
%%%%%%%% END SET UP SPACE AND TIME GRIDS (AND REMAINING PARAMETERS) %%%%%%%%

%%%%%%%% CREATE INITIAL CONDITION %%%%%%%%
% Setup initial conditions
if IC_type==1
    amp_factor_by_k=amp_factor*IC_amps_by_k(wave_number);
    phase_shift=0;
    CA_NL_prepare_disclinsol_for_NLsolver;  % CHECK DEFAULT VALUE FOR NUM_DAMPED!!
    u=umatini*amp_factor;
    v=vmatini*amp_factor;
    z=zmatini*amp_factor;
    thanom=thetamatini*amp_factor;
    th=thanom+thetatildemat;
    qanom=qmatini*amp_factor; %+qrtildemat;
    q=qanom+qbg;
elseif IC_type==2
    % Import them from linear results
    for wave_number=st_wave_number:st_wave_number+numwaves_in_IC-1
        amp_factor_by_k=amp_factor*IC_amps_by_k(wave_number-st_wave_number+1);
        phase_shift=IC_ps_by_k(wave_number-st_wave_number+1);
        CA_NL_prepare_disclinsol_for_NLsolver;  % CHECK DEFAULT VALUE FOR NUM_DAMPED!!
        u=u+umatini*amp_factor_by_k;
        v=v+vmatini*amp_factor_by_k;
        z=z+zmatini*amp_factor_by_k;
        thanom=thanom+thetamatini*amp_factor_by_k;
        qanom=qanom+qmatini*amp_factor_by_k; % NOT ACTUALLY USED IN THE CODE...
    end
    q=qanom+qbg;
    th=thanom+thetatildemat; %+Lv/cp*qv;
elseif IC_type==3
    % Import them from linear results
    for num_damped=num_damped_vec
    for wave_number=st_wave_number:st_wave_number+numwaves_in_IC-1
        amp_factor_by_k=amp_factor*IC_amps_by_k(wave_number-st_wave_number+1);
        phase_shift=IC_ps_by_k(wave_number-st_wave_number+1);
        CA_NL_prepare_disclinsol_for_NLsolver;  % CHECK DEFAULT VALUE FOR NUM_DAMPED!!
        u=u+umatini*amp_factor_by_k;
        v=v+vmatini*amp_factor_by_k;
        z=z+zmatini*amp_factor_by_k;
        thanom=thanom+thetamatini*amp_factor_by_k;
        qanom=qanom+qmatini*amp_factor_by_k; % NOT ACTUALLY USED IN THE CODE...
    end
    end
    q=qanom+qbg;
    th=thanom+thetatildemat; %+Lv/cp*qv;
elseif IC_type==4
    restart_day=100;
    u=usol(:,:,:,restart_day+1);
    v=vsol(:,:,:,restart_day+1);
    z=zsol(:,:,:,restart_day+1);
    th=thsol(:,:,:,restart_day+1);
    thanom=th-thetatildemat;
    q=qsol(:,:,:,restart_day+1);
    qanom=q-qbg;
elseif IC_type==5
    restart_day=100;
    u=usol;
    v=vsol;
    z=zsol;
    th=thsol;
    thanom=th-thetatildemat;
    q=qsol;
    qanom=q-qbg;
end

% Calculate initial p, w from initial dynamical variables
[wmatini,pmatini]=CA_NL_calc_wp(u,v,g,theta0,dz,scalex,scaley,qbgmat,dafrac,qanom,thanom);
p=pmatini;
w=wmatini;
%%%%%%%% END CREATE INITIAL CONDITION %%%%%%%%

%%%%%%%% SAVE INITIAL CONDITION %%%%%%%%
uini=u;
vini=v;
zini=z;
thanomini=thanom;
qanomini=qanom;
qini=q;
thini=th;
pini=p;
wini=w;

curdir=pwd;
%%%%%%%% END SAVE INITIAL CONDITION %%%%%%%%


