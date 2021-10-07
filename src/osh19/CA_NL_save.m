% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 1/11/2017

% Purpose of script:  Save output from results of running FARE_main. 

% Make sure fname_IC is set correctly in the FARE_init file!!

curdir=pwd;

%cd ('/Users/hrogrosky/Documents/CA_NLmodel/CA_NL_data');
zsineamp='t';
if QBGvecversion==2
    fname = 'uclin_verified';
else
    fname=strcat(['CAu',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);    
end
uclin_verified = usol;
mat_fname = strcat(fname, '.mat'); 
tempname = genvarname(fname);
eval([tempname '= uclin_verified;']);
save(mat_fname, tempname);
clear uclin_verified;

if QBGvecversion==2
    fname = 'vclin_verified';
else
    fname=strcat(['CAv',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
vclin_verified = vsol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= vclin_verified;']);
save(mat_fname, tempname);
clear vclin_verified;

if QBGvecversion==2
    fname = 'thanom_verified';
else
    fname=strcat(['CAthanom',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
thanom_verified = thsolanom;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= thanom_verified;']);
save(mat_fname, tempname);
clear thanom_verified;

if QBGvecversion==2
    fname = 'qanom_verified';
else
    fname=strcat(['CAq',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
qanom_verified = qsolanom;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= qanom_verified;']);
save(mat_fname, tempname);
clear qanom_verified;

if QBGvecversion==2
    fname = 'wtot_verified';
else
    fname=strcat(['CAw',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
wtot_verified = wsol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= wtot_verified;']);
save(mat_fname, tempname);
clear wtot_verified;

if QBGvecversion==2
    fname = 'p_verified';
else
    fname=strcat(['CAp',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
p_verified = psol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= p_verified;']);
save(mat_fname, tempname);
clear p_verified;

if QBGvecversion==2
    fname = 'zeta_verified';
else
    fname=strcat(['CAz',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
zeta_verified = zsol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= zeta_verified;']);
save(mat_fname, tempname);
clear zeta_verified;

if QBGvecversion==2
    fname = 'utot_verified';
else
    fname=strcat(['CAut',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
utot_verified = utotsol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= utot_verified;']);
save(mat_fname, tempname);
clear utot_verified;

if QBGvecversion==2
    fname = 'vtot_verified';
else
    fname=strcat(['CAvt',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
vtot_verified = vtotsol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= vtot_verified;']);
save(mat_fname, tempname);
clear vtot_verified;

if QBGvecversion==2
    fname = 'thtot_verified';
else
    fname=strcat(['CAvt',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
thtot_verified = thsol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= thtot_verified;']);
save(mat_fname, tempname);
clear thtot_verified;

if QBGvecversion==2
    fname = 'qtot_verified';
else
    fname=strcat(['CAvt',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
qtot_verified = qsol;
mat_fname=strcat(fname,'.mat'); 
tempname = genvarname(fname);
eval([tempname '= qtot_verified;']);
save(mat_fname, tempname);
clear qtot_verified;

cd(curdir);
