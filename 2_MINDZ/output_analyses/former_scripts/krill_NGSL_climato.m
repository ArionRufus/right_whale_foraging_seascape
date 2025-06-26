% krill_NGSL_climato
%
% Produce figures for the comparison of FTLE results and obsrved accoustic
% data of krill in the GSL
%
% MapsF 2012

%% Import data

krill_load


% Prepare data for spatial analysis

krill_spatial_comp


%% Figure #1: acoustic sampling transects & dense aggregations

% Define the map projection
m_proj('Lambert Conformal Conic','lon',[-70 -55],'lat',[45.5 52])


% Figure assembled with InkScape later...

% Subplot (a)
figure

m_pcolor(lon,lat,bott)

shading interp

m_grid('linestyle','none','fontsize',16)

set(gca,'fontsize',16)

hold on

I = Year==2005;
%m_plot(Lon(I),Lat(I),'marker','.','color',[.5 .5 .5],'markersize',5,'line','none')
m_plot(Lon(I),Lat(I),'marker','.','color',[0 0 1],'markersize',5,'line','none')

[c,h] = m_contour(lon,lat,bottom,[350 350],'color','k','linewidth',.5);

clabel(c,h,'fontsize',10,'labelspacing',500)

[c,h] = m_contour(lon,lat,bottom,[200 200],'color','k','linewidth',1);

clabel(c,h,'fontsize',10,'labelspacing',500)

[c,h] = m_contour(lon,lat,bottom,[100 100],'color','k','linewidth',.5);

clabel(c,h,'fontsize',10,'labelspacing',500)

m_plot(Lon1(Year1==2005&biny==1),Lat1(Year1==2005&biny==1),'o','color',[0 0 1],'markerfacecolor',[0.95 0.95 0.95],'markersize',6,'linewidth',2)
m_plot(Lon1(Year1==2005&biny==1),Lat1(Year1==2005&biny==1),'o','color',[0.95 0.95 0.95],'markerfacecolor',[0.95 0.95 0.95],'markersize',2)

m_text(-56.7,48.5,'Newfoundland','horizontalalignment','center','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_text(-58.7,51.8,'Strait of Belle-Isle','horizontalalignment','center','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_text(-68.6,49.4,'Lower Estuary','horizontalalignment','center','verticalalignment','middle','rotation',30,'fontsize',12,'backgroundcolor','none')

m_text(-63.5,50.6,'Jacques-Cartier Strait','horizontalalignment','center','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_text(-63,49.5,'Anticosti','horizontalalignment','center','verticalalignment','middle','rotation',-30,'fontsize',10,'backgroundcolor','none')

m_text(-63.1,47.6,'Magdalen','horizontalalignment','center','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_text(-63.1,47.3,'Shallows','horizontalalignment','center','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_text(-65.9,48.8,'Gaspe','horizontalalignment','left','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_text(-59.65,47.1,'Cabot','horizontalalignment','left','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_text(-59.25,46.8,'Strait','horizontalalignment','left','verticalalignment','middle','fontsize',12,'backgroundcolor','none')

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-69,51,'(a)','verticalalignment','middle','fontsize',16)

m_text(-62.5,51.5,'2005','horizontalalignment','center','verticalalignment','middle','fontsize',16)

m_ruler([0 0.1],0)

% Print figure
print('-depsc2','-tiff','-painters','Fig01_NEW_a.eps')

print('-djpeg100','Fig01_NEW_a.jpg')

close


% Subplot (b)
figure

m_pcolor(lon,lat,bott)

shading interp

m_grid('linestyle','none','fontsize',16)

set(gca,'fontsize',16)

hold on

I = Year==2007;
m_plot(Lon(I),Lat(I),'marker','.','color',[.5 .5 .5],'markersize',5,'line','none')

[c,h] = m_contour(lon,lat,bottom,[350 350],'color','k','linewidth',.5);

clabel(c,h,'fontsize',10,'labelspacing',500)

[c,h] = m_contour(lon,lat,bottom,[200 200],'color','k','linewidth',1);

clabel(c,h,'fontsize',10,'labelspacing',500)

[c,h] = m_contour(lon,lat,bottom,[100 100],'color','k','linewidth',.5);

clabel(c,h,'fontsize',10,'labelspacing',500)

m_plot(Lon1(Year1==2007&biny==1),Lat1(Year1==2007&biny==1),'kv','markerfacecolor',[0.95 0.95 0.95],'markersize',6,'linewidth',2)
m_plot(Lon1(Year1==2007&biny==1),Lat1(Year1==2007&biny==1),'v','color',[.95 .95 .95],'markerfacecolor',[0.95 0.95 0.95],'markersize',2)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-69,51,'(b)','verticalalignment','middle','fontsize',16)

m_text(-62.5,51.5,'2007','horizontalalignment','center','verticalalignment','middle','fontsize',16)

% Print figure
print('-depsc2','-tiff','-painters','Fig01_NEW_b.eps')

print('-djpeg100','Fig01_NEW_b.jpg')

close


%% Figure #2: histogram of acoustic-derived biomass values for 2005 & 2007

dumm5 = histc(log10(dummy(dummy>0&Year1==2005)),[-inf -3:0.5:3 inf]); dumm5(end) = [];

dumm7 = histc(log10(dummy(dummy>0&Year1==2007)),[-inf -3:0.5:3 inf]); dumm7(end) = [];


figure

x5 = [-3.35:0.5:3.15];

bar(x5,dumm5,0.33,'edgecolor','k','facecolor','k','linewidth',1.5)

axis([-3.5 3.5 0 3.5e3])

set(gca,'fontsize',14)

title('Acoustic-derived biomass distribution')

xlabel('log_{10}(g m^{-2})')

ylabel('Counts')

hold on

x7 = [-3.15:0.5:3.35];

bar(x7,dumm7,0.33,'edgecolor','k','facecolor',[0.95 0.95 0.95],'linewidth',1.5)

legend('2005','2007','location','NE'); legend('boxoff')

plot([-4 2 2 4],[1e6 1e6 -1e6 -1e6],'k--','linewidth',1.5)

annotation('arrow',[.74 .76],[.69 .69],'headlength',6,'headwidth',6,'headstyle','vback1')

text(2.3,2500,'Dense','fontsize',12); text(2.3,2300,'aggregations','fontsize',12);

%text(3.05,180,'2005','fontsize',12)

% Print figure
print('-depsc2','-tiff','-painters','Fig02_NEW.eps')

print('-djpeg100','Fig02_NEW.jpg')

close


%% FTLE principles

if false
    
% Maps; requires m_map toolbox
m_proj('Lambert Conformal Conic','lon',[-68 -61],'lat',[48 50.5])

figure

subplot 211

torigin = datenum(2005,05,01);

t = find(lltime==datenum(2007,8,20)-torigin);

m_plot(squeeze(llon(:,:,t,t)),squeeze(llat(:,:,t,t)),'k.','markersize',3)

m_grid('linestyle','none','fontsize',12)

set(gca,'fontsize',14)

hold on

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1)

m_plot(llon(79,206,t,t),llat(79,206,t,t),'ko','markerfacecolor','b','markersize',6)

m_plot(llon(81,206,t,t),llat(81,206,t,t),'ko','markerfacecolor','y','markersize',6)

m_plot(llon(80,205,t,t),llat(80,205,t,t),'ko','markerfacecolor','c','markersize',6)

m_plot(llon(80,207,t,t),llat(80,207,t,t),'ko','markerfacecolor','g','markersize',6)

m_plot(llon(80,206,t,t),llat(80,206,t,t),'ko','markerfacecolor','r','markersize',6)

title('FTLE_{Backward}  20-Aug-2006','fontsize',16)


subplot 212

m_plot(squeeze(llon(:,:,t,t+10)),squeeze(llat(:,:,t,t+10)),'k.','markersize',3)

m_grid('linestyle','none','fontsize',12)

set(gca,'fontsize',14)

hold on

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1)

m_plot([llon(80,206,t,t+10) llon(80,205,t,t+10)],[llat(80,206,t,t+10) llat(80,205,t,t+10)],'c','linewidth',2)

m_plot(llon(80,205,t,t+10),llat(80,205,t,t+10),'ko','markerfacecolor','c','markersize',6)

m_plot([llon(80,206,t,t+10) llon(80,207,t,t+10)],[llat(80,206,t,t+10) llat(80,207,t,t+10)],'g','linewidth',2)

m_plot(llon(80,207,t,t+10),llat(80,207,t,t+10),'ko','markerfacecolor','g','markersize',6)

m_plot([llon(80,206,t,t+10) llon(79,206,t,t+10)],[llat(80,206,t,t+10) llat(79,206,t,t+10)],'b','linewidth',2)

m_plot(llon(79,206,t,t+10),llat(79,206,t,t+10),'ko','markerfacecolor','b','markersize',6)

m_plot([llon(80,206,t,t+10) llon(81,206,t,t+10)],[llat(80,206,t,t+10) llat(81,206,t,t+10)],'y','linewidth',2)

m_plot(llon(81,206,t,t+10),llat(81,206,t,t+10),'ko','markerfacecolor','y','markersize',6)

m_plot(llon(80,206,t,t+10),llat(80,206,t,t+10),'ko','markerfacecolor','r','markersize',6)

title('FTLE_{Backward}  31-Jul-2006','fontsize',16)


% Print figure
print('-depsc2','-tiff','-painters','Fig_FTLE_principle.eps')

print('-djpeg100','Fig_FTLE_principle.jpg')

close

end

%% Figure #3: example of daily convergence zones

figure

% Maps; requires m_map toolbox
m_proj('Lambert Conformal Conic','lon',[-68 -61],'lat',[48 50.5])

torigin = datenum(2005,05,01);

t = find(lltime==datenum(2007,8,14)-torigin);

FTLEplot = squeeze(Tr(:,:,22)/max(max(Tr(:,:,22)))); % FTLE the 20th of August 2007

m_pcolor(lon,lat,FTLEplot)

set(gca,'fontsize',16)

colormap(gray_inv); caxis([0.2 0.8])

shading flat

m_grid('linestyle','none')

hold on

[c,h] = m_contour(lon,lat,bottom,[350 350],'color','k','linewidth',0.5);

clabel(c,h,'fontsize',10,'labelspacing',500)

[c,h] = m_contour(lon,lat,bottom,[200 200],'color','k','linewidth',2);

clabel(c,h,'fontsize',10,'labelspacing',500)

[c,h] = m_contour(lon,lat,bottom,[100 100],'color','k','linewidth',1);

clabel(c,h,'fontsize',10,'labelspacing',500)


i1 = 78; j1 = 208;

m_line(squeeze(llon(i1,j1,t,[t:t+10])),squeeze(llat(i1,j1,t,[t:t+10])),'color','r','linewidth',1.5)

m_line(squeeze(llon(i1+1,j1,t,[t:t+10])),squeeze(llat(i1+1,j1,t,[t:t+10])),'color','g','linewidth',1.5)

m_line(squeeze(llon(i1-1,j1,t,[t:t+10])),squeeze(llat(i1-1,j1,t,[t:t+10])),'color','b','linewidth',1.5)

m_line(squeeze(llon(i1,j1+1,t,[t:t+10])),squeeze(llat(i1,j1+1,t,[t:t+10])),'color','y','linewidth',1.5)

m_line(squeeze(llon(i1,j1-1,t,[t:t+10])),squeeze(llat(i1,j1-1,t,[t:t+10])),'color','c','linewidth',1.5)

m_plot(llon(i1+1,j1,t,t),llat(i1+1,j1,t,t),'o','markerfacecolor','g','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1-1,j1,t,t),llat(i1-1,j1,t,t),'o','markerfacecolor','b','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1,j1+1,t,t),llat(i1,j1+1,t,t),'o','markerfacecolor','y','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1,j1-1,t,t),llat(i1,j1-1,t,t),'o','markerfacecolor','c','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1,j1,t,t),llat(i1,j1,t,t),'o','markerfacecolor','r','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1,j1,t,t+10),llat(i1,j1,t,t+10),'o','markerfacecolor','r','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1+1,j1,t,t+10),llat(i1+1,j1,t,t+10),'o','markerfacecolor','g','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1-1,j1,t,t+10),llat(i1-1,j1,t,t+10),'o','markerfacecolor','b','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1,j1+1,t,t+10),llat(i1,j1+1,t,t+10),'o','markerfacecolor','y','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1,j1-1,t,t+10),llat(i1,j1-1,t,t+10),'o','markerfacecolor','c','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i1+1,j1,t,t),llat(i1+1,j1,t,t),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1-1,j1,t,t),llat(i1-1,j1,t,t),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1,j1+1,t,t),llat(i1,j1+1,t,t),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1,j1-1,t,t),llat(i1,j1-1,t,t),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1,j1,t,t),llat(i1,j1,t,t),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1,j1,t,t+10),llat(i1,j1,t,t+10),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1+1,j1,t,t+10),llat(i1+1,j1,t,t+10),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1-1,j1,t,t+10),llat(i1-1,j1,t,t+10),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1,j1+1,t,t+10),llat(i1,j1+1,t,t+10),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)

m_plot(llon(i1,j1-1,t,t+10),llat(i1,j1-1,t,t+10),'o','markerfacecolor',[0.95 0.95 0.95],'markeredgecolor',[0.95 0.95 0.95],'markersize',3,'linewidth',1)


i2 = 72; j2 = 202;

m_line(squeeze(llon(i2,j2,t,[t:t+10])),squeeze(llat(i2,j2,t,[t:t+10])),'color','r','linewidth',1.5)

m_line(squeeze(llon(i2+1,j2,t,[t:t+10])),squeeze(llat(i2+1,j2,t,[t:t+10])),'color','g','linewidth',1.5)

m_line(squeeze(llon(i2-1,j2,t,[t:t+10])),squeeze(llat(i2-1,j2,t,[t:t+10])),'color','b','linewidth',1.5)

m_line(squeeze(llon(i2,j2+1,t,[t:t+10])),squeeze(llat(i2,j2+1,t,[t:t+10])),'color','y','linewidth',1.5)

m_line(squeeze(llon(i2,j2-1,t,[t:t+10])),squeeze(llat(i2,j2-1,t,[t:t+10])),'color','c','linewidth',1.5)

m_plot(llon(i2,j2,t,t+10),llat(i2,j2,t,t+10),'o','markerfacecolor','r','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2+1,j2,t,t+10),llat(i2+1,j2,t,t+10),'o','markerfacecolor','g','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2-1,j2,t,t+10),llat(i2-1,j2,t,t+10),'o','markerfacecolor','b','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2,j2+1,t,t+10),llat(i2,j2+1,t,t+10),'o','markerfacecolor','y','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2,j2-1,t,t+10),llat(i2,j2-1,t,t+10),'o','markerfacecolor','c','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2+1,j2,t,t),llat(i2+1,j2,t,t),'o','markerfacecolor','g','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2-1,j2,t,t),llat(i2-1,j2,t,t),'o','markerfacecolor','b','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2,j2+1,t,t),llat(i2,j2+1,t,t),'o','markerfacecolor','y','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2,j2-1,t,t),llat(i2,j2-1,t,t),'o','markerfacecolor','c','markeredgecolor','k','markersize',6,'linewidth',1)

m_plot(llon(i2,j2,t,t),llat(i2,j2,t,t),'o','markerfacecolor','r','markeredgecolor','k','markersize',6,'linewidth',1)


m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

title('FTLEs, {\itT. raschii } DVM behaviour, 14-Aug-2007')


% Print figure
print('-depsc2','-tiff','-painters','Fig03_NEW.eps')

print('-djpeg100','Fig03_NEW.jpg')

close


%% Figure #4: climatology of August 20XX convergence zones

figure

% Maps; requires m_map toolbox
m_proj('Lambert Conformal Conic','lon',[-70 -56],'lat',[45.5 52])

[c1,~] = contour(lon,lat,bottom,[100 100]);

c1(:,c1(1,:)==100) = nan;

[c2,~] = contour(lon,lat,bottom,[200 200]);

c2(:,c2(1,:)==200) = nan;

[c3,~] = contour(lon,lat,bottom,[350 350]);

c3(:,c3(1,:)==350) = nan;

m_pcolor(lon,lat,ftle1)

set(gca,'fontsize',16)

colormap(gray_inv); caxis([0.3 .6])

shading flat

m_grid('linestyle','none')

hold on

m_line(c1(1,:),c1(2,:),'color','k','linewidth',1)
%m_line(c1(1,:),c1(2,:),'color',[.99 .99 .99],'linestyle',':','linewidth',1)

m_line(c2(1,:),c2(2,:),'color','k','linewidth',1.5)
%m_line(c2(1,:),c2(2,:),'color',[.99 .99 .99],'linestyle',':','linewidth',1.5)

m_line(c3(1,:),c3(2,:),'color','k','linewidth',0.5)

m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),mean(U(1:4:end,1:4:end,:),3),mean(V(1:4:end,1:4:end,:),3),4,'color',[.95 .95 .95],'linewidth',1)
m_quiver(lon(1:4:end,1:4:end),lat(1:4:end,1:4:end),mean(U(1:4:end,1:4:end,:),3),mean(V(1:4:end,1:4:end,:),3),4,'color','k','linewidth',0.5)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

title('FTLEs, {\itT. raschii } DVM behaviour, mean of August 2006 to 2011')


% Print figure
print('-depsc2','-tiff','-painters','Fig04_NEW.eps')

print('-djpeg100','Fig04_NEW.jpg')

close


%% Figure #5: observed densities inside zones of high FTLEs at the scale of the GSL

figure

subplot 211

xplot = 0.45+0.05*[1:9]; med = mean(dummy(dummy>0));

plot(xplot,dz(1:end-1)/med,'k-.','linewidth',1.5)

set(gca,'fontsize',14,'xtick',[.5 .55 .6 .65 .7 .75 .8 .85 .9 ],'xticklabel',{'50%','','40%','','30%','','20%','','10%'},'ytick',[1 1.2 1.4 1.6 ])

axis([0.45 0.95 0.8 1.8])

hold on

plot(xplot,ds(1:end-1)/med,'k','linewidth',1.5)

plot(xplot,dx(1:end-1)/med,'k--','linewidth',2)

plot([0 1],[1 1],'k','linewidth',1)

ylabel('Biomass, inside/average')

text(0.9,1.7,'(a)','fontsize',16)

title('Gulf of St. Lawrence','fontsize',16)

legend('Slope','FTLE surface','FTLE DVM','location','northwest')
legend('boxoff')


subplot 212

load NEMO_Gulf_Slope_distance; PZ = double(resnum)./size(result,1);
load Surf_2006-2011_8_BWD_NEW_M2_FTLE_20d_distance; PS = double(resnum)./size(result,1);
load Tr_2006-2011_8_BWD_NEW_M2_FTLE_20d_distance; PX = double(resnum)./size(result,1);

plot(xplot,PZ(1:end-1),'k-.','linewidth',1.5)

set(gca,'fontsize',14,'xtick',[.5 .55 .6 .65 .7 .75 .8 .85 .9],'xticklabel',{'50%','','40%','','30%','','20%','','10%'},'ytick',[.2 .4 .6 .8 ],'yticklabel',{'20%','40%','60%','80%'})

axis([0.45 0.95 0 1])

hold on

plot(xplot,PS(1:end-1),'k','linewidth',2)

plot(xplot,PX(1:end-1),'k--','linewidth',2)

plot([0.45 1],[0.55 0],'k')

ax = [0.595 0.595]; ay = [0.11 0.15];

annotation('arrow',ax,ay,'headlength',6,'headwidth',6,'headstyle','vback1')

xlabel('Proportion of highest values (slope or FTLE) for area selection')

ylabel('Dense aggregations inside')

text(0.9,0.9,'(b)','fontsize',16)


% Print figure
print('-depsc2','-tiff','-painters','Fig05_NEW.eps')

print('-djpeg100','Fig05_NEW.jpg')

close


%% Figure #6: compare distributions of aggregations and slope/FTLE at the scale of the GSL

figure

set(gcf,'papertype','usletter','paperposition',[0 0 8.5 11])

xoff = 0.1; yoff = 0.05;


subplot('position',[0+0.5*xoff 2/3+yoff 1-xoff 1/3-2*yoff])

m_proj('Lambert Conformal Conic','lon',[-70 -56],'lat',[45.5 52])
%m_proj('Lambert Conformal Conic','lon',[-71 -55],'lat',[41 52])

ki = 6; % 25% of highest values for area selection

%distlon = [gulf_dist(:,2,ki) gulf_dist(:,5,ki) nan(size(gulf_dist,1),1)];
%distlon = reshape(distlon',numel(distlon),1);

%distlat = [gulf_dist(:,3,ki) gulf_dist(:,6,ki) nan(size(gulf_dist,1),1)];
%distlat = reshape(distlat',numel(distlat),1);

slopeplot = slope; slopeplot(slope<qz(ki)) = 0; slopeplot(slope>=qz(ki)) = 1;

slopeplot(~isfinite(ftle1(:,:,1))) = nan; 

[c100,~] = contour(lon,lat,bottom,[100 100]); c100(:,c100(1,:)==100) = nan; 

[c200,~] = contour(lon,lat,bottom,[200 200]); c200(:,c200(1,:)==200) = nan;

[c350,~] = contour(lon,lat,bottom,[350 350]); c350(:,c350(1,:)==350) = nan;

m_pcolor(lon,lat,slopeplot)

set(gca,'fontsize',14)

colormap(grey)

shading flat

m_grid('linestyle','none')

hold on

m_line(c100(1,:),c100(2,:),'color','k','linestyle','--','linewidth',1)

m_line(c200(1,:),c200(2,:),'color','k','linestyle','-','linewidth',1.5)

m_line(c350(1,:),c350(2,:),'color','k','linestyle','-','linewidth',.5)

%m_line(distlon,distlat,'color','r','linestyle','-','linewidth',1)

LON=gulf_dist(:,2); LAT=gulf_dist(:,3); YEAR=gulf_dist(:,1); DIST=gulf_dist(:,4);

I = YEAR==2005;
m_plot(LON(I),LAT(I),'ko','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007;
m_plot(LON(I),LAT(I),'kv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2005 & DIST<=4;
m_plot(LON(I),LAT(I),'ro','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007 & DIST<=4;
m_plot(LON(I),LAT(I),'rv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-69,51,'(a)','verticalalignment','middle','fontsize',16)

title(['Inside 25% steepest slope = ' num2str(PZ(ki)*100, '%3.1f') '%' ],'fontsize',16)


subplot('position',[0+0.5*xoff 1/3+yoff 1-xoff 1/3-2*yoff])

%m_proj('Lambert Conformal Conic','lon',[-70 -55],'lat',[45.5 52])

FTLEplot = ftle2; FTLEplot(ftle2<q2(ki)) = 0; FTLEplot(ftle2>=q2(ki)) = 1;

m_pcolor(lon,lat,FTLEplot)

set(gca,'fontsize',14)

colormap(grey)

shading flat

m_grid('linestyle','none')

hold on

m_line(c100(1,:),c100(2,:),'color','k','linestyle','--','linewidth',1)

m_line(c200(1,:),c200(2,:),'color','k','linestyle','-','linewidth',1.5)

m_line(c350(1,:),c350(2,:),'color','k','linestyle','-','linewidth',.5)

LON=Surf_dist(:,2); LAT=Surf_dist(:,3); YEAR=Surf_dist(:,1); DIST=Surf_dist(:,4);

I = YEAR==2005;
m_plot(LON(I),LAT(I),'ko','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007;
m_plot(LON(I),LAT(I),'kv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2005 & DIST<=4;
m_plot(LON(I),LAT(I),'ro','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007 & DIST<=4;
m_plot(LON(I),LAT(I),'rv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-69,51,'(b)','verticalalignment','middle','fontsize',16)

title(['Inside 25% highest FTLEs (surface) = ' num2str(PS(ki)*100, '%3.1f') '%' ],'fontsize',16)


subplot('position',[0+0.5*xoff 0+yoff 1-xoff 1/3-2*yoff])

%m_proj('Lambert Conformal Conic','lon',[-70 -55],'lat',[45.5 52])

dumm1 = ftle1; dumm1(ftle1<q1(ki)) = 0; dumm1(ftle1>=q1(ki)) = 1;

FTLEplot = dumm1;

FTLEplot(FTLEplot>0) = 1;

m_pcolor(lon,lat,FTLEplot)

set(gca,'fontsize',14)

colormap(grey)

shading flat

m_grid('linestyle','none')

hold on

m_line(c100(1,:),c100(2,:),'color','k','linestyle','--','linewidth',1)

m_line(c200(1,:),c200(2,:),'color','k','linestyle','-','linewidth',1.5)

m_line(c350(1,:),c350(2,:),'color','k','linestyle','-','linewidth',.5)

LON=Tr_dist(:,2); LAT=Tr_dist(:,3); YEAR=Tr_dist(:,1); DIST=Tr_dist(:,4);

I = YEAR==2005;
m_plot(LON(I),LAT(I),'ko','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007;
m_plot(LON(I),LAT(I),'kv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2005 & DIST<=4;
m_plot(LON(I),LAT(I),'ro','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007 & DIST<=4;
m_plot(LON(I),LAT(I),'rv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-69,51,'(c)','verticalalignment','middle','fontsize',16)

title(['Inside 25% highest FTLEs (DVM) = ' num2str(PX(ki)*100, '%3.1f') '%' ],'fontsize',16)

% Print figure
print('-depsc2','-tiff','-painters','Fig06_NEW.eps')

print('-djpeg100','Fig06_NEW.jpg')

close


%% Observed densities inside convergence zones around Anticosti Island
%  Require a specific analysis for FTLEs

disp(' ')
disp('Detailed analysis around Anticosti Island')
disp(' ')

% Load data
krill_load

% Compute spatial comparison
krill_spatial_comp

%% Figure #7: observed densities inside zones of high FTLEs around Anticosti

figure

subplot 211

xplot = 0.45+0.05*[1:9]; med = mean(dummy(dummy>0));

plot(xplot,dz(1:end-1)/med,'k-.','linewidth',1.5)

set(gca,'fontsize',14,'xtick',[.5 .55 .6 .65 .7 .75 .8 .85 .9],'xticklabel',{'50%','','40%','','30%','','20%','','10%'},'ytick',[1 1.2 1.4 1.6])

axis([0.45 0.95 0.8 1.8])

hold on

plot(xplot,ds(1:end-1)/med,'k','linewidth',2)

plot(xplot,dx(1:end-1)/med,'k--','linewidth',2)

plot([0 1],[1 1],'k','linewidth',1)

ylabel('Biomass, inside/average')

text(0.9,1.7,'(a)','fontsize',16)

title('Anticosti Island','fontsize',16)

legend('Slope','FTLE surface','FTLE DVM','location','northwest')
legend('boxoff')


subplot 212

load NEMO_Anticosti_Slope_distance; PZ = double(resnum)./size(result,1);
load Surf_2006-2011_8_BWD_NEW_M2_FTLE_Anticosti_20d_distance; PS = double(resnum)./size(result,1);
load Tr_2006-2011_8_BWD_NEW_M2_FTLE_Anticosti_20d_distance; PX = double(resnum)./size(result,1);

plot(xplot,PZ(1:end-1),'k-.','linewidth',1.5)

set(gca,'fontsize',14,'xtick',[.5 .55 .6 .65 .7 .75 .8 .85 .9],'xticklabel',{'50%','','40%','','30%','','20%','','10%'},'ytick',[.2 .4 .6 .8 ],'yticklabel',{'20%','40%','60%','80%'})

axis([0.45 0.95 0 1])

hold on

plot(xplot,PS(1:end-1),'k','linewidth',2)

plot(xplot,PX(1:end-1),'k--','linewidth',2)

plot([0.45 1],[0.55 0],'k')

ax = [0.595 0.595]; ay = [0.11 0.14];

annotation('arrow',ax,ay,'headlength',6,'headwidth',6,'headstyle','vback1')

xlabel('Proportion of highest values (slope or FTLE) for area selection')

ylabel('Dense aggregations inside')

text(0.9,0.9,'(b)','fontsize',16)


% Print figure
print('-depsc2','-tiff','-painters','Fig07_NEW.eps')

print('-djpeg100','Fig07_NEW.jpg')

close


%% Figure #8: compare distributions of aggregations and slope/FTLE around Anticosti

figure

set(gcf,'papertype','usletter','paperposition',[0 0 8.5 11])

xoff = 0.1; yoff = 0.05;


subplot('position',[0+0.5*xoff 2/3+yoff 1-xoff 1/3-2*yoff])

m_proj('Lambert Conformal Conic','lon',[-67 -59],'lat',[48 50.5])

ki = 6; % Percentile for convergence zone = 75%

slopeplot = slope; slopeplot(slope<qz(ki)) = 0; slopeplot(slope>=qz(ki)) = 1;

slopeplot(~isfinite(ftle1(:,:,1))) = nan; 

if ~exist('c100','var')
    [c100,~] = contour(lon,lat,bottom,[100 100]); c100(:,c100(1,:)==100) = nan; 

    [c200,~] = contour(lon,lat,bottom,[200 200]); c200(:,c200(1,:)==200) = nan;

    [c350,~] = contour(lon,lat,bottom,[350 350]); c350(:,c350(1,:)==350) = nan;
end

m_pcolor(lon,lat,slopeplot)

set(gca,'fontsize',14)

colormap(grey)

shading flat

m_grid('linestyle','none')

hold on

m_line([lon(75,213) lon(145,203) lon(145,170) lon(133,170) lon(75,203) lon(75,213)],[lat(75,213) lat(145,203) lat(145,170) lat(133,170) lat(75,203) lat(75,213)],'color','k','linewidth',1)

m_line(c100(1,:),c100(2,:),'color','k','linestyle','--','linewidth',1)

m_line(c200(1,:),c200(2,:),'color','k','linestyle','-','linewidth',1.5)

m_line(c350(1,:),c350(2,:),'color','k','linestyle','-','linewidth',.5)

LON=anti_dist(:,2); LAT=anti_dist(:,3); YEAR=anti_dist(:,1); DIST=anti_dist(:,4);

I = YEAR==2005;
m_plot(LON(I),LAT(I),'ko','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007;
m_plot(LON(I),LAT(I),'kv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2005 & DIST<=4;
m_plot(LON(I),LAT(I),'ro','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007 & DIST<=4;
m_plot(LON(I),LAT(I),'rv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-66,48.5,'(a)','horizontalalignment','center','verticalalignment','bottom','fontsize',16,'backgroundcolor','none')

title(['Inside 25% steepest slope = ' num2str(PZ(ki)*100, '%3.1f') '%' ],'fontsize',16)


subplot('position',[0+0.5*xoff 1/3+yoff 1-xoff 1/3-2*yoff])

%m_proj('Lambert Conformal Conic','lon',[-67 -59],'lat',[48 50.5])

distlon = [Surf_dist(:,1,ki) Surf_dist(:,4,ki) nan(size(Surf_dist,1),1)];
distlon = reshape(distlon',numel(distlon),1);

distlat = [Surf_dist(:,2,ki) Surf_dist(:,5,ki) nan(size(Surf_dist,1),1)];
distlat = reshape(distlat',numel(distlat),1);

FTLEplot = ftle2; FTLEplot(ftle2<q2(ki)) = 0; FTLEplot(ftle2>=q2(ki)) = 1;

m_pcolor(lon,lat,FTLEplot)

set(gca,'fontsize',14)

colormap(grey)

shading flat

m_grid('linestyle','none')

hold on

m_line([lon(75,213) lon(145,203) lon(145,170) lon(133,170) lon(75,203) lon(75,213)],[lat(75,213) lat(145,203) lat(145,170) lat(133,170) lat(75,203) lat(75,213)],'color','k','linewidth',1)

m_line(c100(1,:),c100(2,:),'color','k','linestyle','--','linewidth',1)

m_line(c200(1,:),c200(2,:),'color','k','linestyle','-','linewidth',1.5)

m_line(c350(1,:),c350(2,:),'color','k','linestyle','-','linewidth',.5)

LON=Surf_dist(:,2); LAT=Surf_dist(:,3); YEAR=Surf_dist(:,1); DIST=Surf_dist(:,4);

I = YEAR==2005;
m_plot(LON(I),LAT(I),'ko','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007;
m_plot(LON(I),LAT(I),'kv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2005 & DIST<=4;
m_plot(LON(I),LAT(I),'ro','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007 & DIST<=4;
m_plot(LON(I),LAT(I),'rv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-66,48.5,'(b)','horizontalalignment','center','verticalalignment','bottom','fontsize',16,'backgroundcolor','none')

title(['Inside 25% highest FTLEs (surface) = ' num2str(PS(ki)*100, '%3.1f') '%' ],'fontsize',16)


subplot('position',[0+0.5*xoff 0+yoff 1-xoff 1/3-2*yoff])

%m_proj('Lambert Conformal Conic','lon',[-67 -59],'lat',[48 50.5])

distlon = [Tr_dist(:,1,ki) Tr_dist(:,4,ki) nan(size(Tr_dist,1),1)];
distlon = reshape(distlon',numel(distlon),1);

distlat = [Tr_dist(:,2,ki) Tr_dist(:,5,ki) nan(size(Tr_dist,1),1)];
distlat = reshape(distlat',numel(distlat),1);

dumm1 = ftle1; dumm1(ftle1<q1(ki)) = 0; dumm1(ftle1>=q1(ki)) = 1;

FTLEplot = dumm1;

FTLEplot(FTLEplot>0) = 1;

m_pcolor(lon,lat,FTLEplot)

set(gca,'fontsize',14)

colormap(grey)

shading flat

m_grid('linestyle','none')

hold on

m_line([lon(75,213) lon(145,203) lon(145,170) lon(133,170) lon(75,203) lon(75,213)],[lat(75,213) lat(145,203) lat(145,170) lat(133,170) lat(75,203) lat(75,213)],'color','k','linewidth',1)

m_line(c100(1,:),c100(2,:),'color','k','linestyle','--','linewidth',1)

m_line(c200(1,:),c200(2,:),'color','k','linestyle','-','linewidth',1.5)

m_line(c350(1,:),c350(2,:),'color','k','linestyle','-','linewidth',.5)

LON=Tr_dist(:,2); LAT=Tr_dist(:,3); YEAR=Tr_dist(:,1); DIST=Tr_dist(:,4);

I = YEAR==2005;
m_plot(LON(I),LAT(I),'ko','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007;
m_plot(LON(I),LAT(I),'kv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2005 & DIST<=4;
m_plot(LON(I),LAT(I),'ro','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'o','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

I = YEAR==2007 & DIST<=4;
m_plot(LON(I),LAT(I),'rv','markerfacecolor',[.95 .95 .95],'markersize',6,'linewidth',1.5)
m_plot(LON(I),LAT(I),'v','color',[.95 .95 .95],'markerfacecolor',[.95 .95 .95],'markersize',3)

m_line(NEMO_lon,NEMO_lat,'color','k','linewidth',1.5)

m_text(-66,48.5,'(c)','horizontalalignment','center','verticalalignment','bottom','fontsize',16,'backgroundcolor','none')

title(['Inside 25% highest FTLEs (DVM) = ' num2str(PX(ki)*100, '%3.1f') '%' ],'fontsize',16)


% Print figure
print('-depsc2','-tiff','-painters','Fig08_NEW.eps')

print('-djpeg100','Fig08_NEW.jpg')

close
