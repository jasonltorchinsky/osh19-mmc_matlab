% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 2/6/2017
% Script name:  make_timelong_figs

% Purpose of script:  make x-t Hovmoller plots of results of nonlinear FARE
% model 
 
vartoplot=1; % 1-u, 2-v, 3-p, 4-w, 5-th, 6-q
tnumtoskip=1;
levtoplot=3; % 1-bottom of troposphere, numlevs - top of troposphere

% if vartoplot<=3 levtoplot=numlevs+1;
% else 
%     if mod(numlevs,2)==1
%         levtoplot=(numlevs+1)/2;
%     else
%         levtoplot=(numlevs+2)/2;
%     end
% end

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

numstoplot=1:tnumtoskip:tnumtosave;
filenamestringbase1=strcat(['CAt']); 
if QBGvecversion==2
    filenamestringbase2=strcat(['yEQz',num2str(round(zzU(levtoplot)*10)/10),'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
else
    filenamestringbase2=strcat(['yEQz',num2str(round(zzU(levtoplot)*10)/10),'Bvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
end

xxrs=[xx*360/40000 360];
[XXtl TTtl]=meshgrid(xxrs,tts(numstoplot)/3600/24);

figure;
set(gcf,'Position',[200,200,400,550]);
hold on;
% colormap(jet);
colormap(redblue(21));
qmax=max(max(abs((squeeze(qsol(numy/2,:,levtoplot,numstoplot(end)))-squeeze(qbg(numy/2,:,levtoplot)))))); %transpose on squeeze(qsol?

contourf(XXtl,TTtl,[transpose(squeeze(mean(qsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1)))-repmat(squeeze(mean(qbg(numy/2:numy/2+1,:,levtoplot),1)),[length(numstoplot),1]) squeeze(mean(qsol(numy/2:numy/2+1,1,levtoplot,numstoplot),1))-repmat(squeeze(mean(qsol(numy/2:numy/2+1,1,levtoplot),1)),[length(numstoplot),1])],[-qmax:qmax/30:qmax],'Linestyle','none');

colorbar('Location','EastOutside');
caxis([-qmax qmax]);
maxth=max(max(abs(transpose(squeeze(mean(thsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1)))-repmat(squeeze(mean(thetatildemat(numy/2:numy/2+1,:,levtoplot),1)),[length(numstoplot),1]))));
% maxth=max(max(max(abs(squeeze(thsolanom(numy/2,:,:,numstoplot(end))))))); %-thetatildemat(numy/2,:,:))))));
thinc=maxth/3;

contour(XXtl,TTtl,[transpose(squeeze(mean(thsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1)))-repmat(squeeze(mean(thetatildemat(numy/2:numy/2+1,:,levtoplot),1)),[length(numstoplot),1]) (squeeze(mean(thsol(numy/2:numy/2+1,1,levtoplot,numstoplot),1)))-repmat(squeeze(mean(thetatildemat(numy/2:numy/2+1,1,levtoplot))),[length(numstoplot),1])],[thinc/2:thinc:maxth-thinc/2],'Linestyle','-','Color','k'); %-squeeze(mean(thetatildemat(numy/2:numy/2+1,:,levtoplot),1))   % -thetatildemat(:,1,levtoplot)]
contour(XXtl,TTtl,[transpose(squeeze(mean(thsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1)))-repmat(squeeze(mean(thetatildemat(numy/2:numy/2+1,:,levtoplot),1)),[length(numstoplot),1]) (squeeze(mean(thsol(numy/2:numy/2+1,1,levtoplot,numstoplot),1)))-repmat(squeeze(mean(thetatildemat(numy/2:numy/2+1,1,levtoplot))),[length(numstoplot),1])],[-(maxth-thinc/2):thinc:-thinc/2],'Linestyle','--','Color','k'); %-squeeze(mean(thetatildemat(numy/2:numy/2+1,:,levtoplot),1))   % -thetatildemat(:,1,levtoplot)]

% if vartoplot==1
%     contour(XXtl,TTtl,1000*[transpose(squeeze(mean(utotsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1))) squeeze(mean(utotsol(numy/2:numy/2+1,1,levtoplot,numstoplot),1))],10);
%     title(strcat(['u (m/s), EQ, z=',num2str(round(zzU(levtoplot)*10)/10),' km']));
%     filenamestring=strcat([filenamestringbase1,'u',fnameIC,filenamestringbase2]);
% elseif vartoplot==2
%     contour(XXtl,TTtl,1000*[transpose(squeeze(mean(vtotsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1))) squeeze(mean(vtotsol(numy/2:numy/2+1,1,levtoplot,numstoplot),1))],10);
%     title(strcat(['v (m/s), EQ, z=',num2str(round(zzU(levtoplot)*10)/10),' km']));
%     filenamestring=strcat([filenamestringbase1,'v',fnameIC,filenamestringbase2]);
% elseif vartoplot==3
%     contour(XXrs,ZZrs,[transpose(squeeze(mean(psol(numy/2:numy/2+1,:,2:end,numtoplot),1))) squeeze(mean(psol(numy/2:numy/2+1,1,2:end,numtoplot),1))],10);
%     title(strcat(['p, EQ']));
%     filenamestring=strcat([filenamestringbase1,'p',fnameIC,filenamestringbase2]);
% elseif vartoplot==4
%     contour(XXrs,ZZrs,1000*[transpose(squeeze(mean(wsol(numy/2:numy/2+1,:,:,numtoplot),1))) squeeze(mean(wsol(numy/2:numy/2+1,1,:,numtoplot),1))],10);
%     title(strcat(['w (m/s), EQ']));
%     filenamestring=strcat([filenamestringbase1,'w',fnameIC,filenamestringbase2]);
% elseif vartoplot==5
%     contour(XXtl,TTtl,[transpose(squeeze(mean(thsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1))) squeeze(mean(thsol(numy/2:numy/2+1,1,levtoplot,numstoplot),1))],10);
%     title(strcat(['\theta, EQ, z=',num2str(round(zzW(levtoplot)*10)/10),' km']));
%     filenamestring=strcat([filenamestringbase1,'th',fnameIC,filenamestringbase2]);
% elseif vartoplot==6
%     contour(XXtl,TTtl,[transpose(squeeze(mean(qsol(numy/2:numy/2+1,:,levtoplot,numstoplot),1))) squeeze(mean(qsol(numy/2:numy/2+1,1,levtoplot,numstoplot),1))],10);
%     title(strcat(['q, EQ, z=',num2str(round(zzW(levtoplot)*10)/10),' km']));
%     filenamestring=strcat([filenamestringbase1,'q',fnameIC,filenamestringbase2]);
% end
xlim([0 360]);
ylim([0 tts(numstoplot(end))/3600/24]);
% ylim([0 10]);
yl=ylim;
ylimmax=yl(2);
set(gca,'XTick',[0 90 180 270 360]);
set(gca,'XTickLabel',{'0','90','180','270','360'});
% set(gca,'YTick',[0:4:16]);
ylabel('t (days)');
xlabel('Longitude (deg)');
% text(10,15,strcat(['t=',num2str(round(tts(numtoplot)/3600/24*10)/10)]));
colorbar('Location','EastOutside');

hold on;
num_damped=1;
speedtoplot=Frequencies(num_damped)/(kbar_int/40000000);
% angletorotate=70;
offset=28;
numlines=1;
for pcount=0:numlines-1
%     plot([0 360]+360/numlines*pcount+offset,[0 40000000/speedtoplot],'--k','Linewidth',1);
end
angletorotate=55;
stpms=40000000/(40000000/speedtoplot*86400);
% text(2+offset,tts(numstoplot(end))/3600/24*0.04,strcat([num2str(round(10*stpms)/10),' m/s']),'rotation',angletorotate);

hold on;
speedtoplot=Frequencies(num_damped)/(kbar_int/40000000)*1.4;
% angletorotate=70;
offset=0;
% toffset=150;
% numlines=1;
% for pcount=0:numlines-1
%     plot([0 360]+360/numlines*pcount+offset,[toffset toffset+40000000/speedtoplot],'-k','Linewidth',1);
% end
plot([0,0],[0 200],'k');
plot([360,360],[0 200],'k');
plot([0,360],[0 0],'k');
plot([0,360],[200 200],'k');
toffset=110;
for pcount=0:numlines-1
%     plot([0 360]+360/numlines*pcount+offset,[toffset toffset+40000000/speedtoplot],'-k','Linewidth',1);
end
angletorotate=44;
stpms=40000000/(40000000/speedtoplot*86400);
% text(10,68+toffset-59,strcat([num2str(round(10*stpms)/10),' m/s']),'rotation',angletorotate,'Fontweight','Bold');

figureHandle = gcf;
set(gca,'fontsize',figfontsize);
set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
set(gcf,'Units','points');
set(gcf,'PaperUnits','points');
sizepaper=get(gcf,'Position');
sizepaper=sizepaper(3:4);
set(gcf,'PaperSize',sizepaper);
set(gcf,'PaperPosition',[0,0,sizepaper(1),sizepaper(2)]);

filenamestring=strcat([filenamestringbase1,fnameIC,filenamestringbase2]);
filenamestring=strrep(filenamestring,'.','p');
filenamestring=strrep(filenamestring,'-','n');
% saveas(gcf,filenamestring,'epsc');

cd(curdir);
