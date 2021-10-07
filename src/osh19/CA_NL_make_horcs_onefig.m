% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 2/6/2017
% Script name:  make_horcs_figs

% Purpose of script:  plot x-y snapshots of results of nonlinear FARE model 

% close all;

% Load data 
vartoplot=1;   % 1-u, 2-v, 3-p, 4-w, 5-th, 6-q - NOTE, 3 AND 4 MAY NOT WORK
plotevery=20;  % Make a plot at t=0, plotevery, 2*plotevery, etc.
levtoplot=7;

% % Set level that we want to see.  
% % DEFAULT = highest level for u,v,p; level near middle of the troposphere for theta_e,q_t,w
% if vartoplot<=3 levtoplot=3; %3; %10; %numlevs+1; 
% else
%     if mod(numlevs,2)==1
%         levtoplot=(numlevs-1)/2;
%     else
%         levtoplot=(numlevs)/2;
%     end
% end

% Filename for saving figure
filenamestringbase1=strcat(['CAh']);
if QBGvecversion==2
    filenamestringbase2=strcat(['C7z',num2str(round(zzU(levtoplot)*10)/10),'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
else    
    filenamestringbase2=strcat(['C7z',num2str(round(zzU(levtoplot)*10)/10),'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end
% curdir=pwd;
% 
% cd ('/Users/hrogrosky/Documents/CA_NLmodel/CA_NL_data');
% 
% if QBGvecversion==2
%     fname=strcat(['CAu',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAu',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['usol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAv',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAv',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['vsol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAth',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAth',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['thsol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAq',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAq',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['qsol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAw',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAw',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['wsol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAp',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAp',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['psol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAz',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAz',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['zsol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAut',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAut',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['utotsol=' tempname ';']);
% 
% if QBGvecversion==2
%     fname=strcat(['CAvt',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% else
%     fname=strcat(['CAvt',fnameIC,'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'te',num2str(tend),'dtp',num2str(dtstart*saveevery),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
% end
% fname=strrep(fname,'.','p');
% fname=strrep(fname,'-','n');
% mat_fname=strcat(fname,'.mat'); 
% load(mat_fname);
% tempname = genvarname(fname);
% eval(['vtotsol=' tempname ';']);
% 
% cd(curdir);
% 
% cd ('/Users/hrogrosky/Documents/VCU_papers/CA_NL_report/CA_figs');

XXrs=[XX*360/40000 XX(:,1)+360];
YYrs=[YY*360/40000 YY(:,1)*360/40000];
nxskip=2;
nyskip=2;
XXrssp=XXrs(1:nyskip:end,1:nxskip:end);
YYrssp=YYrs(1:nyskip:end,1:nxskip:end);

numtoplotrange=1:plotevery:tnumtosave; %40:161; %191; %241; %41; %1, 41, 81, 121, 161 %161:plotevery:201;

for numtoplot=numtoplotrange %161:plotevery:241 %216 %301 %201 %113 %tend/3600/24+1   % Could also set the last day for plotting manually.  %101 %5:51 %5:51 %1:11 %51
    Usparse=1000*[(squeeze(utotsol(1:nyskip:end,1:nxskip:end,levtoplot,numtoplot))) squeeze(utotsol(1:nyskip:end,1,levtoplot,numtoplot))]; %transpose(squeeze... on first term?
    Vsparse=1000*[(squeeze(vtotsol(1:nyskip:end,1:nxskip:end,levtoplot,numtoplot))) squeeze(vtotsol(1:nyskip:end,1,levtoplot,numtoplot))]; %transpose(squeeze... on first term?

    figure;
    set(gcf,'Position',[200,200,500,200]); %[200,200,500,180] or ... 500 250]

    quivscale=4;
    hold on;
%     colormap(jet);
    colormap(redblue(21));
    qmax=max(max(abs((squeeze(max(qsol(:,:,levtoplot,numtoplotrange),[],4)))-squeeze(qbg(:,:,levtoplot))))); %transpose on squeeze(qsol?

    contourf(XXrs,YYrs,[qsol(:,:,levtoplot,numtoplot)-qbg(:,:,levtoplot) qsol(:,1,levtoplot,numtoplot)-qbg(:,1,levtoplot)],[-qmax:qmax/30:qmax],'Linestyle','none'); % 
    colorbar('Location','EastOutside');
    caxis([-qmax qmax]);
    
    quiver(XXrssp,YYrssp,Usparse/quivscale,Vsparse/quivscale,'k','Linewidth',1.25);
    maxth=max(max(abs(squeeze(max(thsol(:,:,levtoplot,numtoplotrange),[],4)-thetatildemat(:,:,levtoplot)))));
%     maxth=max(max(max(abs(squeeze(thsol(:,:,:,numtoplot)-thetatildemat)))));
    thinc=maxth/5;
    
    hold on;
    contour(XXrs,YYrs,[thsol(:,:,levtoplot,numtoplot)-thetatildemat(:,:,levtoplot) thsol(:,1,levtoplot,numtoplot)-thetatildemat(:,1,levtoplot)],[thinc/2:thinc:maxth-thinc/2],'Linestyle','-','Color','k'); 
    contour(XXrs,YYrs,[thsol(:,:,levtoplot,numtoplot)-thetatildemat(:,:,levtoplot) thsol(:,1,levtoplot,numtoplot)-thetatildemat(:,1,levtoplot)],[-(maxth-thinc/2):thinc:-thinc/2],'Linestyle','--','Color','k'); 
    
    plot([0 0],[-45 45],'k');
    plot([360 360],[-45 45],'k');
    plot([0 360],[-45 -45],'k');
    plot([0 360],[45 45],'k');
    
    filenamestring=strcat([filenamestringbase1,fnameIC,'t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
%     if vartoplot==1
%         maxu=max(max(abs(1000*utotsol(:,:,levtoplot,numtoplot))));
%         maxtomin=2*maxu/19;
%         contour(XXrs,YYrs,1000*[squeeze(utotsol(:,:,levtoplot,numtoplot)) squeeze(utotsol(:,1,levtoplot,numtoplot))],[-maxu:maxtomin:maxu]); %,[-3:0.5:-1 1:0.5:3]);
%         title(strcat(['u_{tot} (m/s), z=',num2str(round(zzU(levtoplot)*10)/10),' km']));
%         filenamestring=strcat([filenamestringbase1,'ut',fnameIC,'t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% %         caxis([-3 3]);
%     elseif vartoplot==2
%         contour(XXrs,YYrs,1000*[squeeze(vtotsol(:,:,levtoplot,numtoplot)) squeeze(vtotsol(:,1,levtoplot,numtoplot))]); %,[-4:-1 1:4]);
%         title(strcat(['v_{tot} (m/s), z=',num2str(round(zzU(levtoplot)*10)/10),' km']));
%         filenamestring=strcat([filenamestringbase1,'vt',fnameIC,'t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% %         caxis([-4 4]);
%     elseif vartoplot==3
%         contour(XXrs,YYrs,squeeze(psol(:,:,levtoplot,numtoplot)),10);
%         title(strcat(['p, z=',num2str(zzU(levtoplot)),' km']));
%         filenamestring=strcat([filenamestringbase1,'p',fnameIC,'t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
%     elseif vartoplot==4
%         contour(XXrs,YYrs,1000*squeeze(wsol(:,:,levtoplot,numtoplot)),10);
%         title(strcat(['w (m/s), z=',num2str(zzW(levtoplot)),' km']));
%         filenamestring=strcat([filenamestringbase1,'w',fnameIC,'t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
%     elseif vartoplot==5
%         contour(XXrs,YYrs,[squeeze(thsol(:,:,levtoplot,numtoplot))-thetatildemat(:,:,levtoplot) squeeze(thsol(:,1,levtoplot,numtoplot))-thetatildemat(:,1,levtoplot)]); %,[-1:0.1:-0.1 0.1:0.1:1]);
%         title(strcat(['\theta anom (K), z=',num2str(round(zzW(levtoplot)*10)/10),' km']));
%         filenamestring=strcat([filenamestringbase1,'th',fnameIC,'t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% %         caxis([-1 1]);
%     elseif vartoplot==6
%         contour(XXrs,YYrs,[squeeze(qsol(:,:,levtoplot,numtoplot))-qbg(:,:,levtoplot) squeeze(qsol(:,1,levtoplot,numtoplot))-qbg(:,1,levtoplot)]); %,[-0.003:0.0005:-0.0005 0.0005:0.0005:0.003]);
%         title(strcat(['q anom (kg/kg), z=',num2str(round(zzW(levtoplot)*10)/10),' km']));
%         filenamestring=strcat([filenamestringbase1,'q',fnameIC,'t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% %         caxis([-0.003 0.003]);
%     end
    xlim([0 360]);
%     ylim([-pY*360/40000,pY*360/40000]);
    ylim([-45 45]);
    set(gca,'XTick',[0 90 180 270 360]);
    set(gca,'XTickLabel',{'0','90','180','270','360'});
    set(gca,'YTick',[-40 -20 0 20 40]);
    set(gca,'YTickLabel',{'40S','20S','EQ','20N','40N'});
    
%     text(10,pY*360/40000*0.8,strcat(['t=',num2str(round(tts(numtoplot)/3600/24*10)/10)]));
    text(10,45*0.8,strcat(['t=',num2str(round(tts(numtoplot)/3600/24*10)/10)]));
    colorbar('Location','EastOutside');
    zstring=num2str(round(10*zzU(levtoplot))/10); 
%     title(strcat(['z=',zstring,' km']));
    xlabel('Longitude (deg)');
    if numtoplot==1 text(-40,-60,'(a)');
    elseif numtoplot==41 text(-40,-60,'(c)');
    elseif numtoplot==81 text(-40,-60,'(e)');
    elseif numtoplot==121 text(-40,-60,'(g)');
%     elseif numtoplot==161 text(-40,-60,'(i)');
    end
    if numtoplot==161 text(-40,-60,'(a)'); end
    figureHandle = gcf;
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    set(gcf,'Units','points');
    set(gcf,'PaperUnits','points');
    sizepaper=get(gcf,'Position');
    sizepaper=sizepaper(3:4);
    set(gcf,'PaperSize',sizepaper);
    set(gcf,'PaperPosition',[0,0,sizepaper(1),sizepaper(2)]);

    filenamestring=strrep(filenamestring,'.','p');
    filenamestring=strrep(filenamestring,'-','n');
%     saveas(gcf,filenamestring,'epsc');

    horstrucwindamp=sqrt(Usparse.^2+Vsparse.^2);
% disp('---------------------------------------------------------------- ');
% disp(strcat(['t=',int2str(numtoplot-1)]));
% disp('Maximum wind in horizontal structure plot is ');
% max(max(abs(horstrucwindamp)))
% disp('Maximum theta anomaly magnitude in horizontal structure plot is ');
% maxth
% disp('theta contour interval drawn at the following fractions of maxth ');
% [thinc/2:thinc:maxth-thinc/2]
% [-(maxth-thinc/2):thinc:-thinc/2]
% disp('---------------------------------------------------------------- ');

end
% 
% horstrucwindamp=sqrt(Usparse.^2+Vsparse.^2);
% disp('---------------------------------------------------------------- ');
% disp(strcat(['t=',int2str(numtoplot-1)]));
% disp('Maximum wind in horizontal structure plot is ');
% max(max(abs(horstrucwindamp)))
% disp('Maximum theta anomaly magnitude in horizontal structure plot is ');
% maxth
% disp('theta contour interval drawn at the following fractions of maxth ');
% [thinc/2:thinc:maxth-thinc/2]
% [-(maxth-thinc/2):thinc:-thinc/2]
% disp('---------------------------------------------------------------- ');

cd(curdir);
