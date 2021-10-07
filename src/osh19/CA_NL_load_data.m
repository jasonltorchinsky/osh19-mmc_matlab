% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 4/27/2017

% Purpose of script:  Use to load a previously-computed IC

curdir=pwd;

%%%%%%%% ONLY USE THE LINES HERE IF YOU NEED TO RE-USE AN INITIAL CONDITION
%%%%%%%% FROM A DIFFERENT SET OF PARAMETER VALUES %%%%%%%%
% thdamp_tau_temp=thdamp_tau;
% thdamp_tau=30*3600*24;
% VTtemp=VTdim_mps;
% VTdim_mps=0.02;
% Bvsdimtemp=Bvsdim;
% Bvsdim=-0.00128;
% numlevstemp=numlevs;
% numlevs=40;
% Dhu_temp=Dhu;
% Dht_temp=Dht;
% Dhq_temp=Dhq;
% Dhu=0.3;
% Dht=0.3;
% Dhq=0.3;
% numxtemp=numx;
% numx=40;
% numytemp=numy;
% numy=40;

%%%%%%%% END (BUT ALSO SEE BELOW!!) ONLY USE THE LINES HERE IF YOU NEED TO
%%%%%%%% RE-USE AN INITIAL CONDITION FROM A DIFFERENT SET OF PARAMETER 
%%%%%%%% VALUES %%%%%%%%

cd ('/Users/hrogrosky/Documents/FARE_NL_rotation_BTBC/FARE_data_BTBC');

fname=strcat(['FNLui_Bvs',num2str(Bvsdim),'_VT',num2str(VTdim_mps),'_Y',int2str(pY),'_x',int2str(numx),'_y',int2str(numy),'_z',int2str(numlevs),'_Dhu',num2str(Dhu),'_Dht',num2str(Dht),'_Dhq',num2str(Dhq),'_tauth',num2str(thdamp_tau)]);
fname=strrep(fname,'.','p');
fname=strrep(fname,'-','n');
mat_fname=strcat(fname,'.mat'); 
load(mat_fname);
tempname = genvarname(fname);
eval(['usol=' tempname ';']);
fnamevars={fname};
clear(fnamevars{:});

fname=strcat(['FNLvi_Bvs',num2str(Bvsdim),'_VT',num2str(VTdim_mps),'_Y',int2str(pY),'_x',int2str(numx),'_y',int2str(numy),'_z',int2str(numlevs),'_Dhu',num2str(Dhu),'_Dht',num2str(Dht),'_Dhq',num2str(Dhq),'_tauth',num2str(thdamp_tau)]);
fname=strrep(fname,'.','p');
fname=strrep(fname,'-','n');
mat_fname=strcat(fname,'.mat'); 
load(mat_fname);
tempname = genvarname(fname);
eval(['vsol=' tempname ';']);
fnamevars={fname};
clear(fnamevars{:});

fname=strcat(['FNLtei_Bvs',num2str(Bvsdim),'_VT',num2str(VTdim_mps),'_Y',int2str(pY),'_x',int2str(numx),'_y',int2str(numy),'_z',int2str(numlevs),'_Dhu',num2str(Dhu),'_Dht',num2str(Dht),'_Dhq',num2str(Dhq),'_tauth',num2str(thdamp_tau)]);
fname=strrep(fname,'.','p');
fname=strrep(fname,'-','n');
mat_fname=strcat(fname,'.mat'); 
load(mat_fname);
tempname = genvarname(fname);
eval(['tesol=' tempname ';']);
fnamevars={fname};
clear(fnamevars{:});

fname=strcat(['FNLqti_Bvs',num2str(Bvsdim),'_VT',num2str(VTdim_mps),'_Y',int2str(pY),'_x',int2str(numx),'_y',int2str(numy),'_z',int2str(numlevs),'_Dhu',num2str(Dhu),'_Dht',num2str(Dht),'_Dhq',num2str(Dhq),'_tauth',num2str(thdamp_tau)]);
fname=strrep(fname,'.','p');
fname=strrep(fname,'-','n');
mat_fname=strcat(fname,'.mat'); 
load(mat_fname);
tempname = genvarname(fname);
eval(['qtsol=' tempname ';']);
fnamevars={fname};
clear(fnamevars{:});

fname=strcat(['FNLwi_Bvs',num2str(Bvsdim),'_VT',num2str(VTdim_mps),'_Y',int2str(pY),'_x',int2str(numx),'_y',int2str(numy),'_z',int2str(numlevs),'_Dhu',num2str(Dhu),'_Dht',num2str(Dht),'_Dhq',num2str(Dhq),'_tauth',num2str(thdamp_tau)]);
fname=strrep(fname,'.','p');
fname=strrep(fname,'-','n');
mat_fname=strcat(fname,'.mat'); 
load(mat_fname);
tempname = genvarname(fname);
eval(['wsol=' tempname ';']);
fnamevars={fname};
clear(fnamevars{:});

fname=strcat(['FNLpi_Bvs',num2str(Bvsdim),'_VT',num2str(VTdim_mps),'_Y',int2str(pY),'_x',int2str(numx),'_y',int2str(numy),'_z',int2str(numlevs),'_Dhu',num2str(Dhu),'_Dht',num2str(Dht),'_Dhq',num2str(Dhq),'_tauth',num2str(thdamp_tau)]);
fname=strrep(fname,'.','p');
fname=strrep(fname,'-','n');
mat_fname=strcat(fname,'.mat'); 
load(mat_fname);
tempname = genvarname(fname);
eval(['psol=' tempname ';']);
fnamevars={fname};
clear(fnamevars{:});

fname=strcat(['FNLzi_Bvs',num2str(Bvsdim),'_VT',num2str(VTdim_mps),'_Y',int2str(pY),'_x',int2str(numx),'_y',int2str(numy),'_z',int2str(numlevs),'_Dhu',num2str(Dhu),'_Dht',num2str(Dht),'_Dhq',num2str(Dhq),'_tauth',num2str(thdamp_tau)]);
fname=strrep(fname,'.','p');
fname=strrep(fname,'-','n');
mat_fname=strcat(fname,'.mat'); 
load(mat_fname);
tempname = genvarname(fname);
eval(['zsol=' tempname ';']);
fnamevars={fname};
clear(fnamevars{:});

%%%%%%%% ONLY USE THE LINES HERE IF YOU NEED TO RE-USE AN INITIAL CONDITION
%%%%%%%% FROM A DIFFERENT SET OF PARAMETER VALUES %%%%%%%%
thdamp_tau=thdamp_tau_temp;
VTdim_mps=VTtemp;
Bvsdim=Bvsdimtemp;
numlevs=numlevstemp;
Dhu=Dhu_temp;
Dht=Dht_temp;
Dhq=Dhq_temp;
numx=numxtemp;
numy=numytemp;
% fnameIC=fnameICnew;
%%%%%%%% END ONLY USE THE LINES HERE IF YOU NEED TO RE-USE AN INITIAL 
%%%%%%%% CONDITION FROM A DIFFERENT SET OF PARAMETER VALUES %%%%%%%%

cd(curdir);

% Check to see if we need to add/subtract x, y points
usol=FARE_resize_IC(usol,numx,numy,numlevs+1,H,1);
vsol=FARE_resize_IC(vsol,numx,numy,numlevs+1,H,1);
zsol=FARE_resize_IC(zsol,numx,numy,1,H,1);
qtsol=FARE_resize_IC(qtsol,numx,numy,numlevs+1,H,2);
tesol=FARE_resize_IC(tesol,numx,numy,numlevs+1,H,2);
psol=FARE_resize_IC(psol,numx,numy,numlevs+1,H,1);
wsol=FARE_resize_IC(wsol,numx,numy,numlevs+1,H,2);
