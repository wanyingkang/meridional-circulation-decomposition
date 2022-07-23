casenames={'HS_S5','HotPole_ColdEQ_S5','HS_season','HotPoleColdEQ_season'};
cranges=[5e-9 1e-9 1e-8 2e-9];
crange1s=[50 20 100 100];

showtop=50;
showbot=1000;
logscale=0;

%for icase=1:length(casenames)
for icase=2:2
    
    casename=char(casenames(icase))
    crange=cranges(icase);
    crange1=crange1s(icase);
    
% casename='HS_S5';
% crange=4e-9;
% crange1=100;
 
% casename='HotPole_ColdEQ_S5';
% crange=1e-9;
% crange1=10;

% DJF
% casename='HS_season';
% crange=1e-8;
% crange1=100;

%  casename='HotPoleColdEQ_season';
%  crange=2e-9;
% crange1=100;

fileclim=['./All_',casename,'.nc'];
pltalls=[false,false];

%% read
nbot0=0;
lat=ncread(fileclim,'lat');
nlat=length(lat);
lev=ncread(fileclim,'lev');
nlev=length(lev);
VU=squeeze(mean(ncread(fileclim,'VU',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
VT=squeeze(mean(ncread(fileclim,'VT',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
V=squeeze(mean(ncread(fileclim,'V',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
U=squeeze(mean(ncread(fileclim,'U',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
T=squeeze(mean(ncread(fileclim,'T',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
Teq=squeeze(mean(ncread(fileclim,'TREF',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
OMEGA=squeeze(mean(ncread(fileclim,'OMEGA',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
OMEGAU=squeeze(mean(ncread(fileclim,'OMEGAU',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
OMEGAT=squeeze(mean(ncread(fileclim,'OMEGAT',[1,1,1,1],[Inf,Inf,nlev-nbot0,1])));
PS=squeeze(mean(ncread(fileclim,'PS',[1,1,1],[Inf,Inf,1])));
lev=lev(1:end-nbot0);
nlev=nlev-nbot0;

% lev_=repmat(lev'.*100,[nlat,1]);
% PS_=repmat(PS',[1,nlev]);
% mask=(lev_>PS_);
% mask(1,:)=true;mask(end,:)=true;
 VU(1,:)=0.;VU(end,:)=0.;
 VT(1,:)=0.;VT(end,:)=0.;
 V(1,:)=0.;V(end,:)=0.;


%% calculate omega equation terms
vt_eddy=VT-V.*T;
omegat_eddy=OMEGAT-OMEGA.*T;
uv_eddy=VU-V.*U;
omegau_eddy=OMEGAU-U.*OMEGA;

f0=2*2*pi/86400;
g=9.81;
a=6371000;
R=287;
kappa=2./7.;
f=repmat(reshape(f0*sind(lat),[nlat,1]),[1,nlev]);
clat=repmat(reshape(cosd(lat),[nlat,1]),[1,nlev]);
i_lat=[1:nlat];
p_lat=min(i_lat+1,nlat);
m_lat=max(i_lat-1,1);
dlat=repmat(reshape(lat(p_lat)-lat(m_lat),[nlat,1]),[1,nlev]).*(pi/180);
i_lev=[1:nlev];
p_lev=min(i_lev+1,nlev);
m_lev=max(i_lev-1,1);
z=log(1000./lev);
dz=repmat(reshape(z(p_lev)-z(m_lev),[1,nlev]),[nlat,1]);
dlev=repmat(reshape(lev(p_lev)-lev(m_lev),[1,nlev]),[nlat,1]).*100;
Th_T=(1000./lev').^kappa;

tmp=1./(clat.^2).*(uv_eddy(p_lat,:).*clat(p_lat,:).^2 - uv_eddy(m_lat,:).*clat(m_lat,:).^2)./dlat./a;
%tmp=1./(clat).*(uv_eddy(p_lat,:).*clat(p_lat,:) - uv_eddy(m_lat,:).*clat(m_lat,:))./dlat./a;
Term_uv=-f.*(tmp(:,p_lev)-tmp(:,m_lev))./dz;
Term_uv=Term_uv.*clat;

UV_y=1./clat.*V.*(U(p_lat,:).*clat(p_lat,:)-U(m_lat,:).*clat(m_lat,:))./dlat./a;
Term_UVadv=-f.*(UV_y(:,p_lev)-UV_y(:,m_lev))./dz;
Term_UVadv=Term_UVadv.*clat;

tmp=(omegau_eddy(:,p_lev) - omegau_eddy(:,m_lev))./dlev;
Term_uomega=-f.*(tmp(:,p_lev)-tmp(:,m_lev))./dz;
Term_uomega=Term_uomega.*clat;

tmp=OMEGA.*(U(:,p_lev) - U(:,m_lev))./dlev;
Term_WUadv=-f.*(tmp(:,p_lev)-tmp(:,m_lev))./dz;
Term_WUadv=Term_WUadv.*clat;

ks=1/86400/4; ka=1/86400/40;
sigmab=0.7;
k=zeros(size(f));
ibelow=lev>sigmab*1000;
k(:,ibelow)=ka+(ks-ka).*(clat(:,ibelow).^4).*(lev(ibelow)'/1000-sigmab)./(1-sigmab);
k(:,~ibelow)=ka;
J=k.*(Teq-T);
Term_Q=R.*(J(p_lat,:)-J(m_lat,:))./dlat./a;

kf=1/86400;
kuv=zeros(size(f));
kuv(:,ibelow)=kf.*(lev(ibelow)'/1000-sigmab)./(1-sigmab).*ones(nlat,1);
FU=-kuv.*U;
FV=-kuv.*V;
Term_FU=f.*(FU(:,p_lev)-FU(:,m_lev))./dz;
Term_FU=Term_FU.*clat;

tmp=1./clat.*(vt_eddy(p_lat,:).*clat(p_lat,:) - vt_eddy(m_lat,:).*clat(m_lat,:))./dlat./a;
tmp(1,:)=0;tmp(end,:)=0;
Term_vt=-R.*(tmp(p_lat,:)-tmp(m_lat,:))./dlat./a;
Term_vt=Term_vt.*clat;

tmp=V.*(T(p_lat,:) - T(m_lat,:))./dlat./a;
Term_VTadv=-R.*(tmp(p_lat,:)-tmp(m_lat,:))./dlat./a;
Term_VTadv=Term_VTadv.*clat;

tmp=1./Th_T.*(omegat_eddy(:,p_lev).*Th_T(p_lev)-omegat_eddy(:,m_lev).*Th_T(m_lev))./dlev;
Term_omegat=-R.*(tmp(p_lat,:)-tmp(m_lat,:))./dlat./a;

tmp=1./Th_T.*OMEGA.*(T(:,p_lev).*Th_T(p_lev)-T(:,m_lev).*Th_T(m_lev))./dlev;
Term_WTadv=-R.*(tmp(p_lat,:)-tmp(m_lat,:))./dlat./a;

Term_ADV=Term_VTadv+Term_WUadv+Term_UVadv;

Th=T.*(Th_T);
S=-T.*(log(Th(:,p_lev))-log(Th(:,m_lev)))./dlev;
Term_stm=f.^2.*(V(:,p_lev)-V(:,m_lev))./dz + R.*(S(p_lat,:).*OMEGA(p_lat,:)-S(m_lat,:).*OMEGA(m_lat,:))./dlat./a;
Term_stm=-Term_stm.*clat;

% Term_uv(mask)=0.;
% Term_UVadv(mask)=0.;
% Term_uomega(mask)=0.;
% Term_WUadv(mask)=0.;
% Term_Q(mask)=0.;
% Term_FU(mask)=0.;
% Term_vt(mask)=0.;
% Term_VTadv(mask)=0.;
% Term_omegat(mask)=0.;
% Term_WTadv(mask)=0.;
% Term_ADV(mask)=0.;
% Term_stm(mask)=0.;

%% plot omega equation terms
for ip=1:length(pltalls)
    
pltall=pltalls(ip);
if pltall
    totalplt=9;
else
    totalplt=5;
end

figure
polarmap;
iplt=1;
subplot(1,totalplt,iplt)
contourf(lat,lev,Term_uv',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange+crange/20 crange-crange/20])
set(gca,'yscale','log','ydir','reverse')
title('eddy uv')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
ylabel('Pressure (mb)')
xlabel('latitude')
iplt=iplt+1;
subplot(1,totalplt,iplt)
contourf(lat,lev,Term_vt',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('eddy vt')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat,lev,Term_uomega',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('eddy u\omega')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

if pltall
subplot(1,totalplt,iplt)
contourf(lat,lev,Term_omegat',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('\omega t')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;
subplot(1,totalplt,iplt)
contourf(lat,lev,Term_Q',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('Q')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;
subplot(1,totalplt,iplt)
contourf(lat,lev,Term_FU',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('friction')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;
subplot(1,totalplt,iplt)
contourf(lat,lev,Term_ADV',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('Advection')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat,lev,Term_uv'+Term_uomega'+Term_Q'+Term_vt'+Term_omegat'+Term_ADV'+Term_FU',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('Sum')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;
else
subplot(1,totalplt,iplt)
contourf(lat,lev,Term_uv'+Term_uomega'+Term_vt',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('Sum')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;
end

subplot(1,totalplt,iplt)
contourf(lat,lev,Term_stm',[-crange*10,[-crange-crange/20:crange/10:crange+crange/20],crange*10],'LineColor','none')
caxis([-crange-crange/20 crange+crange/20])
set(gca,'yscale','log','ydir','reverse')
title('\nabla^2\Psi_m')
xticks([-90:30:90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

hp4=[0.55    0.08    0.33    0.13];
colorbar('Position', [hp4(1)+hp4(3)+0.04  hp4(2)+ hp4(4)-0.02 0.007  hp4(4)*5],'Ticks',[-crange,-crange/2,0,crange/2,crange])

if pltall
    set(gcf,'Position',[237         633        2500         315])
    saveas(gcf,['./laplace_psi_',casename,'_all.png'])
else
    set(gcf,'Position',[237         633        1800         315])
    saveas(gcf,['./laplace_psi_',casename,'.png'])
end
end
%% solve poisson
if ~isfile(['psi_driver_',casename,'.mat'])
global nlat nlev lat lev clat S f R a
   [psi_uv0,lat1,lev1]=poisson_stm(-Term_uv);
   [psi_vt0,~,~]=poisson_stm(-Term_vt);
   [psi_uomega0,lat1,lev1]=poisson_stm(-Term_uomega);
   [psi_omegat0,~,~]=poisson_stm(-Term_omegat);
   [psi_Q0,~,~]=poisson_stm(-Term_Q);
   [psi_FU0,~,~]=poisson_stm(-Term_FU);
   [psi_ADV0,~,~]=poisson_stm(-Term_ADV);
%% meridional streamfunction
ilev=ncread(fileclim,'ilev');ilev=ilev(1:end-nbot0);
ilev1=ilev(2:end);
dp=100.*(ilev(2:end)-ilev(1:end-1));
psi0=2*pi*a/g*cosd(lat).*cumsum(V.*dp',2);
 
%% unit
unit=1/g*2*pi*a/1e9;
 psi_uv=psi_uv0.*unit;
 psi_vt=psi_vt0.*unit;
 psi_uomega=psi_uomega0.*unit;
 psi_omegat=psi_omegat0.*unit;
 psi_Q=psi_Q0.*unit;
 psi_FU=psi_FU0.*unit;
 psi_ADV=psi_ADV0.*unit;
 
 psi1=psi0/1e9;
 [zz,yy]=meshgrid(ilev1,lat);
 [zz1,yy1]=meshgrid(lev1,lat1);
 psi=interp2(zz,yy,psi1,zz1,yy1);
  save(['psi_driver_',casename,'.mat'],'psi','psi_FU','psi_Q','psi_uomega','psi_omegat','psi_uv','psi_vt','psi_ADV','lev1','lat1');
  %save(['psi_driver_',casename,'.mat'],'psi_ADV','-append')
else
  load(['psi_driver_',casename,'.mat'])
end

psi_FU=smooth2a(psi_FU,1,2);
psi=smooth2a(psi,5,1);
 %% plot
 totalplt=5;
 iplt=1;

figure
polarmap

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uv'+psi_vt'+psi_uomega'+psi_omegat',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('eddy')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_FU',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('friction')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_Q',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('\psi_Q')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uv'+psi_vt'+psi_uomega'+psi_omegat'+psi_Q'+psi_FU',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('Sum')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi'-psi_ADV',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('\Psi-\Psi_{adv}')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

hp4=[0.55    0.08    0.33    0.13];
colorbar('Position', [hp4(1)+hp4(3)+0.04  hp4(2)+ hp4(4)-0.02 0.009  hp4(4)*5.5],'Ticks',[-crange1,-crange1/2,0,crange1/2,crange1])


set(gcf,'Position',[237         633        1800         315])
%saveas(gcf,['./psi_decomp_',casename,'.png'])
saveas(gcf,['~/Desktop/psi_decomp_',casename,'.png'])

%%
figure
polarmap
totalplt=4;
iplt=1;
subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uv',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('\psi_{uv}')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_vt',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('\psi_{vt}')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uomega',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('\psi_{u\omega}')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_FU',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('friction')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

hp4=[0.55    0.08    0.33    0.13];
colorbar('Position', [hp4(1)+hp4(3)+0.04  hp4(2)+ hp4(4)-0.02 0.009  hp4(4)*5.5],'Ticks',[-crange1,-crange1/2,0,crange1/2,crange1])

set(gcf,'Position',[237         633        1800         315])
%saveas(gcf,['./psi_decomp_eddy_',casename,'.png'])
saveas(gcf,['~/Desktop/psi_decomp_eddy_',casename,'.png'])

%%
figure
polarmap
totalplt=5;
iplt=1;

if icase==1 || icase==2

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uv',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==1
title('\psi_{uv}')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uomega',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==1
title('\psi_{u\omega}')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_FU',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==1
title('friction')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_FU'+psi_uv'+psi_uomega',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==1
title('sum')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==1
title('\psi')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;
end

if icase == 3 || icase==4
subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uv'+psi_vt',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==3
title('\psi_{uv}+\psi_{vt}')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_FU',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==3
title('friction')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_Q',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==3
title('diabatic')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi_uv'+psi_vt'+psi_FU'+psi_uv'+psi_Q',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==3
title('sum')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;

subplot(1,totalplt,iplt)
contourf(lat1,lev1,psi',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
if icase ==3
title('\psi')
end
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
iplt=iplt+1;
end
hp4=[0.55    0.08    0.33    0.13];
colorbar('Position', [hp4(1)+hp4(3)+0.04  hp4(2)+ hp4(4)-0.02 0.009  hp4(4)*5.5],'Ticks',[-crange1,-crange1/2,0,crange1/2,crange1])

set(gcf,'Position',[237         633        1800         315])
%saveas(gcf,['./psi_decomp_select_',casename,'.png'])
saveas(gcf,['~/Desktop/psi_decomp_select_',casename,'.png'])


%% ADV
figure
polarmap
contourf(lat1,lev1,psi_ADV',[-crange1*10,[-crange1-crange1/20:crange1/10:crange1+crange1/20],crange1*10],'LineColor','none')
caxis([-crange1-crange1/20 crange1-crange1/20])
set(gca,'yscale','log','ydir','reverse')
title('\psi_{ADV}')
xticks([-90:30:90])
xlim([-90 90])
yticks([4,10,20,50,100,200,500,950])
ylim([showtop,showbot])
set(gca,'FontSize',20)
xlabel('latitude')
end

%% Poisson solver
% function [psi,lat1,lev1]=poisson_stm(T)
% global nlat nlev lat lev clat S f R a mask
% nbot=0;
% lev1=lev(1:end-nbot);
% P0=101325;
% pres=lev1.*100;pres=cat(1,pres,P0);
% lat1=lat;
% clat1=clat(:,1:end-nbot);clat1=cat(2,clat1,clat1(:,end));
% seclat1=1./clat1;seclat1(1,:)=1.;seclat1(end,:)=1.;
% S1=S(:,1:end-nbot);S1=cat(2,S1,S1(:,end));
% f1=f(:,1:end-nbot);f1=cat(2,f1,f1(:,end));
% T=T(:,1:end-nbot);T=cat(2,T,T(:,end));
% mask1=cat(2,mask,mask(:,end));mask1(:,end)=true;
% d2r=pi/180;
% 
% psi0=zeros(nlat,nlev-nbot+1);
% psi=psi0;
% dt0=2e8;
% dt=dt0;
% err=1;
% 
% 
% for it=1:100000
% if err>1e-13
%     psi0=psi;
%     dt=dt0/(it)^0.1;
%     [psi_p,psi_y]=gradient(psi,pres,lat1.*d2r);
%     [psi_pp,~]=gradient(psi_p,pres,lat1.*d2r);
%     [~,psi_yy]=gradient(psi_y.*S1.*seclat1,pres,lat1.*d2r);
%     tend=dt.*(f1.^2.*pres'.*psi_pp+R/a^2.*clat1.*psi_yy+T);
%     %tend(1,:)=0;tend(end,:)=0;
%     tend(:,1)=0;%tend(:,end)=0;
%     tend(mask1)=0;
%     psi=psi0+tend;
%     err=max(tend(:)/dt);
%     disp(it);disp(err)
% end
% end
% psi=psi(:,1:end-1);
% end

% %% Poisson solver
function [psi,lat1,lev1]=poisson_stm(T)
global nlat nlev lat lev clat S f R a
nbot=0;
lev1=lev(1:end-nbot);
pres=lev1.*100;
P0=101325;
pres=lev1.*100;pres=cat(1,pres,P0);
lat1=lat;
clat1=clat(:,1:end-nbot);clat1=cat(2,clat1,clat1(:,end));
seclat1=1./clat1;seclat1(1,:)=seclat1(2,:);seclat1(end,:)=seclat1(end-1,:);
S1=S(:,1:end-nbot);S1=cat(2,S1,S1(:,end));
f1=f(:,1:end-nbot);f1=cat(2,f1,f1(:,end));
T=T(:,1:end-nbot);T=cat(2,T,T(:,end));
d2r=pi/180;

psi0=zeros(nlat,nlev-nbot+1);
psi=psi0;
dt0=8e8;
dt=dt0;
err=1;


for it=1:200000
if err>1e-16
    psi0=psi;
    dt=dt0/it^0.1;
    [psi_p,psi_y]=gradient(psi,pres,lat1.*d2r);
    [psi_pp,~]=gradient(psi_p,pres,lat1.*d2r);
    [~,psi_yy]=gradient(psi_y.*S1.*seclat1,pres,lat1.*d2r);
    tend=dt.*(f1.^2.*pres'.*psi_pp+R/a^2.*clat1.*psi_yy+T);
    tend(1,:)=0;tend(end,:)=0;tend(:,1)=0;tend(:,end)=0;
    psi=psi0+tend;
    err=max(tend(:)./dt);
    if mod(it,1000)==0
    disp(it);disp(err)
    end
end
end

psi=psi(:,1:end-1);
end