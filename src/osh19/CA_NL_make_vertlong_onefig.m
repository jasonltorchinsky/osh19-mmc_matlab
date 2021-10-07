% Reed Ogrosky
% Virginia Commonwealth University, Department of Mathematics and Applied
% Mathematics
% 2/6/2017
% Script name:  make_vertlong_figs

% Purpose of script:  plot x-z snapshots of results of nonlinear FARE model 

% close all;

vartoplot=1;   % 1-u, 2-v, 3-p, 4-w, 5-th, 6-q - NOTE, 3 AND 4 MAY NOT WORK
plotevery=20;  % Make a plot at t=0, plotevery, 2*plotevery, etc.

% curdir=pwd;
% 
% cd ('/Users/hrogrosky/Documents/CA_NLmodel/CA_NL_data');

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

% Make plots

XXrs=[transpose(squeeze(XX3U(1,:,:)))*360/40000 (squeeze(XX3U(1,1,:)))+360];
nxskip=2;
nyskip=2;
XXrssp = XXrs(:,1:nxskip:end);
% ZZu = [dz/2:dz:H-dz/2];    %YYrs(2:nyskip:48,1);
%     ZZw = [dz:dz:H-dz];
% ZZw = [0:dz:H];
% ZZurssp = ZZu;
% ZZwrssp = ZZw;

ZZurs=[transpose(squeeze(ZZ3U(1,:,:))) (squeeze(ZZ3U(1,1,:)))]; 
ZZurssp=ZZurs(2:end,1:nxskip:end); XXrsshort=XXrs(2:end,:); XXrsspshort=XXrssp(2:end,:);
ZZwrs=[transpose(squeeze(ZZ3W(1,:,:))) (squeeze(ZZ3W(1,1,:)))]; 
ZZwrssp=ZZwrs(:,1:nxskip:end);
% ZZurssp=ZZurs;
% ZZwrssp=ZZwrs;
% ZZwrs=ZZurs(2:end,:); 

numtoplotrange=1:plotevery:tnumtosave;

for numtoplot=numtoplotrange %191 %241 %1:plotevery:401 %161:plotevery:201 %tend/3600/24+1 %:plotevery:134 %tend/3600/24+1 %51
    Usparse=1000*[transpose(squeeze(mean(utotsol(numy/2:numy/2+1,1:nxskip:end,2:end,numtoplot),1))) squeeze(mean(utotsol(numy/2:numy/2+1,1,2:end,numtoplot),1))];
    Wsparse=1000*[transpose(squeeze(mean(wsol(numy/2:numy/2+1,1:nxskip:end,2:end,numtoplot),1))) squeeze(mean(wsol(numy/2:numy/2+1,1,2:end,numtoplot),1))];
    
    filenamestringbase1=strcat(['CAv']); 
    if QBGvecversion==2
        filenamestringbase2=strcat(['yEQBvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'exs',num2str(expscale),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
    else
        filenamestringbase2=strcat(['yEQBvs',num2str(Bvsdim),'Y',int2str(pY),'x',int2str(numx),'y',int2str(numy),'z',int2str(numlevs),'dt',num2str(dtstart),'zsa',num2str(zsineamp),'tt',num2str(tautop),'tb',num2str(taubottom),'bt',num2str(btop),'bb',num2str(bbottom),'tu',num2str(tau_u),'tth',num2str(tau_theta),'Dvt',num2str(Dvtop),'Dvb',num2str(Dvbottom)]);
    end
    figure;
    
    set(gcf,'Position',[200,200,500,200]);
    vartoplot=1;
    quivscale=2;
    hold on;
    colormap(jet);
    colormap(redblue(21));
    qmax=max(max(abs(transpose(squeeze(mean(max(qsol(numy/2:numy/2+1,:,:,numtoplotrange),[],4),1))-squeeze(mean(qbg(numy/2:numy/2+1,:,:),1))))));

    contourf(XXrs,ZZwrs,[transpose(squeeze(mean(max(qsol(numy/2:numy/2+1,:,:,numtoplot),[],4),1))-squeeze(mean(qbg(numy/2:numy/2+1,:,:),1))) squeeze(mean(qsol(numy/2:numy/2+1,1,:,numtoplot),1))-squeeze(mean(qbg(numy/2:numy/2+1,1,:),1))],[-qmax:qmax/30:qmax],'Linestyle','none'); %,[-0.003:0.0005:-0.0005 0.0005:0.0005:0.003]);
    colorbar('Location','EastOutside');
    caxis([-qmax qmax]);
%     caxis([-0.000001 0.000001]);
    hq=quiver(XXrsspshort,ZZurssp,Usparse/quivscale,150*Wsparse/quivscale,'k','Maxheadsize',0.01,'Linewidth',1.25); %'ShowArrowHead','off'); %1/quivscale
    
    maxth=max(max(abs(squeeze(min(thsol(numy/2,:,:,numtoplotrange),[],4)-thetatildemat(numy/2,:,:)))));
    thinc=maxth/5;
    
    hold on;
    contour(XXrs,ZZwrs,[transpose(squeeze(mean(thsol(numy/2:numy/2+1,:,:,numtoplot),1))-squeeze(mean(thetatildemat(numy/2:numy/2+1,:,:),1))) squeeze(mean(thsol(numy/2:numy/2+1,1,:,numtoplot),1))-squeeze(mean(thetatildemat(numy/2:numy/2+1,1,:),1))],[thinc/2:thinc:maxth-thinc/2],'Linestyle','-','Color','k'); %,[-1:0.1:-0.1 0.1:0.1:1]);
    
    hold on;
    contour(XXrs,ZZwrs,[transpose(squeeze(mean(thsol(numy/2:numy/2+1,:,:,numtoplot),1))-squeeze(mean(thetatildemat(numy/2:numy/2+1,:,:),1))) squeeze(mean(thsol(numy/2:numy/2+1,1,:,numtoplot),1))-squeeze(mean(thetatildemat(numy/2:numy/2+1,1,:),1))],[-(maxth-thinc/2):thinc:-thinc/2],'Linestyle','--','Color','k'); %,[-1:0.1:-0.1 0.1:0.1:1]);
    
    plot([0,360],[0,0],'k');
    plot([0,360],[16,16],'k');
    plot([360,360],[0,16],'k');
    plot([0,0],[0,16],'k');

%     headWidth = 8;
%     headLength = 3;
%     LineLength = 0.08;
% 
%     U = hq.UData;
%     V = hq.VData;
%     X = hq.XData;
%     Y = hq.YData;
% 
%     [nxt nyt]=size(X);
%     for ii = 1:nxt
%         for ij = 1:nyt
% 
%             headWidth = 2;
%             ah = annotation('arrow',...
%                 'HeadLength',headLength,'HeadWidth',headWidth);
%             set(ah,'parent',gca);
%             set(ah,'position',[X(ii,ij) Y(ii,ij) LineLength*U(ii,ij) LineLength*V(ii,ij)]);
% 
%         end
%     end
    hold on;
    
%     contour(XXrs,ZZrs,1000*[transpose(squeeze(mean(utotsol(numy/2:numy/2+1,:,2:end,numtoplot),1))) squeeze(mean(utotsol(numy/2:numy/2+1,1,2:end,numtoplot),1))]); %,[-10:-1 1:10]);
%         title(strcat(['u (m/s), EQ']));
%         filenamestring=strcat([filenamestringbase1,'u',fnameIC,'_t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
%         caxis([-10 10]);
% % %     elseif vartoplot==2
% % %         contour(XXrs,ZZrs,1000*[transpose(squeeze(mean(vtotsol(numy/2+4:numy/2+1+4,:,2:end,numtoplot),1))) squeeze(mean(vtotsol(numy/2:numy/2+1,1,2:end,numtoplot),1))]); %,[-10:-1 1:10]);
% % %         title(strcat(['v (m/s), ~10N']));
% % %         filenamestring=strcat([filenamestringbase1,'v',fnameIC,'_t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% % % %         caxis([-10 10]);
% % %     elseif vartoplot==3
% % %         contour(XXrs,ZZrs,[transpose(squeeze(mean(psol(numy/2:numy/2+1,:,2:end,numtoplot),1))) squeeze(mean(psol(numy/2:numy/2+1,1,2:end,numtoplot),1))],10);
% % %         title(strcat(['p, EQ']));
% % %         filenamestring=strcat([filenamestringbase1,'p',fnameIC,'_t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% % %     elseif vartoplot==4
% % %         contour(XXrs,ZZrs,1000*[transpose(squeeze(mean(wsol(numy/2:numy/2+1,:,:,numtoplot),1))) squeeze(mean(wsol(numy/2:numy/2+1,1,:,numtoplot),1))],10);
% % %         title(strcat(['w (m/s), EQ']));
% % %         filenamestring=strcat([filenamestringbase1,'w',fnameIC,'_t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% % %     elseif vartoplot==5
% % %         contour(XXrs,ZZrs,[transpose(squeeze(mean(thsol(numy/2:numy/2+1,:,:,numtoplot),1))-squeeze(mean(thetatildemat(numy/2:numy/2+1,:,:),1))) squeeze(mean(thsol(numy/2:numy/2+1,1,:,numtoplot),1))-squeeze(mean(thetatildemat(numy/2:numy/2+1,1,:),1))]); %,[-1:0.1:-0.1 0.1:0.1:1]);
% % %         title(strcat(['\theta anom (K), EQ']));
% % %         filenamestring=strcat([filenamestringbase1,'th',fnameIC,'_t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% % % %         caxis([-1 1]);
% % %     elseif vartoplot==6
% % %         contour(XXrs,ZZrs,[transpose(squeeze(mean(qsol(numy/2:numy/2+1,:,:,numtoplot),1))-squeeze(mean(qbg(numy/2:numy/2+1,:,:),1))) squeeze(mean(qsol(numy/2:numy/2+1,1,:,numtoplot),1))-squeeze(mean(qbg(numy/2:numy/2+1,1,:),1))]); %,[-0.003:0.0005:-0.0005 0.0005:0.0005:0.003]);
% % %         title(strcat(['q anom (kg/kg), EQ']));
% % %         filenamestring=strcat([filenamestringbase1,'q',fnameIC,'_t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
% % % %         caxis([-0.003 0.003]);
% % %     end
    xlim([0 360]);
    ylim([0 H]);
    set(gca,'XTick',[0 90 180 270 360]);
    set(gca,'XTickLabel',{'0','90','180','270','360'});
    set(gca,'YTick',[0:4:16]);
    ylabel('z (km)');
    if numtoplot==1 text(-40,-1,'(b)');
    elseif numtoplot==41 text(-40,-1,'(d)');
    elseif numtoplot==81 text(-40,-1,'(f)');
    elseif numtoplot==121 text(-40,-1,'(h)');
%     elseif numtoplot==161 text(-40,-1,'(j)');
    end
    if numtoplot==161 text(-40,-1,'(b)'); end
    text(10,14.5,strcat(['t=',num2str(round(tts(numtoplot)/3600/24*10)/10)]));
%     title(strcat(['y=',ystring]));
    xlabel('Longitude (deg)');
%     colorbar('Location','EastOutside');
    
    figureHandle = gcf;
    set(gca,'fontsize',figfontsize);
    set(findall(figureHandle,'type','text'),'fontSize',figfontsize)
    set(gcf,'Units','points');
    set(gcf,'PaperUnits','points');
    sizepaper=get(gcf,'Position');
    sizepaper=sizepaper(3:4);
    set(gcf,'PaperSize',sizepaper);
    set(gcf,'PaperPosition',[0,0,sizepaper(1),sizepaper(2)]);
    
    cd ('/Users/hrogrosky/Documents/VCU_papers/CA_NL_report/CA_figs');
    filenamestring=strcat([filenamestringbase1,fnameIC,'_t',num2str(round(tts(numtoplot)/3600/24*10)/10),filenamestringbase2]);
    
    filenamestring=strrep(filenamestring,'.','p');
    filenamestring=strrep(filenamestring,'-','n');
    saveas(gcf,filenamestring,'epsc');

%     vertstrucwindamp=sqrt(Usparse.^2+Wsparse.^2);
% disp('---------------------------------------------------------------- ');
% disp(strcat(['t=',int2str(numtoplot-1)]));
% disp('Maximum wind in vertical structure plot is ');
% max(max(abs(vertstrucwindamp)))
% disp('Maximum theta anomaly magnitude in vertical structure plot is ');
% maxth
% disp('theta contour interval drawn at the following fractions of maxth ');
% [thinc/2:thinc:maxth-thinc/2]
% [-(maxth-thinc/2):thinc:-thinc/2]
% disp('---------------------------------------------------------------- ');

end

% vertstrucwindamp=sqrt(Usparse.^2+Wsparse.^2);
% disp('---------------------------------------------------------------- ');
% disp(strcat(['t=',int2str(numtoplot-1)]));
% disp('Maximum wind in vertical structure plot is ');
% max(max(abs(vertstrucwindamp)))
% disp('Maximum theta anomaly magnitude in vertical structure plot is ');
% maxth
% disp('theta contour interval drawn at the following fractions of maxth ');
% [thinc/2:thinc:maxth-thinc/2]
% [-(maxth-thinc/2):thinc:-thinc/2]
% disp('---------------------------------------------------------------- ');

% cd(curdir);
