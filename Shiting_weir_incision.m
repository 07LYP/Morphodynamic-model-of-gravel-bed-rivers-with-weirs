clear;
clc;
tic;
%load the pre-earthquake equilibrium surface texture
load('PreEarthquake_WC.mat');

%%%%%%Input of Model Parameters%%%%%%
%ep: small number;
%If: flood intermittency;
%timeyear: coefficient convert year to second;
%g: gravitaional coefficient;
%rho: water density;
%R: submerged specific gravity;
%L: channel length;
%SI: initial slope;
%Duration: duration of calculation (year);
%qwf: flow discharge per unit width at the entrance (m2/s);
%qtfT: total sediment supply per unit width (m2/s);
%alpr: coefficient in M-S resistance formula;
%nk: coefficient relating roughness height;
%Phi: phi value of grain size;
%D: characteristic grain size of each size range;
%ng: number of sediment group;
%qgfT: gravel supply rate per unit width (m2/s);
%qsfT: sand supply rate per unit width (m2/s);
%qtfT: total sediment supply rate per unit width (m2/s);
%qtf: sediment supply rate per unit width for each range (m2/s);
%fgravel: grain size distribution of gravel supply;
%fsand: grain size distribution of sand supply;
%fsupply: grain size distribution of total supply;
%fbase: grain size distribution of substrate base;
%na: coefficient relating active layer thickness;
%p: bed porosity;
%au: coefficient related to Exner discretization;
%Ls: thickness of storage layer;
%sigma: subsidence rate (m/year);
ep=1e-16;
If=0.03;
timeyear=If*365.25*24*3600;
g=9.8;
rho=1000;
R=1.65;
L=50000;
Duration=20;
Nt=200;
qwf=6;
alpr=8.1;
alpha=0.5;
nk=2;
Phin=linspace(-2,7,10);
Phi=(Phin(1:9)+Phin(2:10))/2;
D=2.^Phi/1000;
ng=size(D,2);
fgravel=[0,0,0,10,10,15,30,25,10]/100;
fsand=[20,30,50,0,0,0,0,0,0]/100;
qgfT=2e-3;
qsfT=2e-3;
qtf=qgfT*fgravel+qsfT*fsand;
qtfT=sum(qtf,2);
fsupply=qtf./qtfT;
fbase=fgravel;
na=2;
p=0.35;
au=1;
Ls=1;

%%%%%%Spatial grid length and temporal step%%%%%%
%dt: time step (year);
%dx: cell size;
%n: cell number;
dt=0.0005;
dx=100;
n=L/dx+1;
n_dam = 5000/dx+1;                  % grid number of the weir
h_dam = 5;                          % weir height
Courant=zeros(n,1);
interval=round(Duration/dt/Nt);

%%%%%%Initial Conditions%%%%%%
%Fi: proportion of sediment on bed surface;
%FS: proportion of sand on bed surface;
%Dsg: geometrical mean grain size of bed surface;
%D90:grain size that 90 percent is finer;
%La: thickness of active layer;
Fi = ones(n,1)*fbase1;
FS = sum(Fi(:,1:3),2);
Dsg= 2.^(Fi*Phi')/1000;
D90=ones(n,1);
Cdg=[zeros(n,1),cumsum(Fi,2)];
for iD=1:n
    ind=find(Cdg(iD,:)>=0.9,1);
    D1=Phin(ind-1);
    D2=Phin(ind);
    P1=Cdg(iD,ind-1);
    P2=Cdg(iD,ind);
    D90(iD,1)=2^((0.9-P1)/(P2-P1)*(D2-D1)+D1)/1000;
end
La  = na*D90;
La0 = La(1);
etaa = 0;       % the alluvial thickness at the weir grid
%zb: bed elevation;
%zbini: initial bed elevation;
%s: bed slope;
x =linspace(0,L,n)';
zb=interp1([0, L], [zb(1), 0], x);
zb(n_dam) = zb(n_dam)+h_dam;
zbini=zb;
s=[(zb(1:n-1)-zb(2:n))/dx;(zb(n-1)-zb(n))/dx];
t=0;
i=0;
%qw: flow discharge per unit width;
%h: water depth;
qw  = qwf*ones(n,1);
h   = zeros(n,1);
% Backwater Equation
dx1 = 1;
n1  = L/dx1+1;
x1  = linspace(0,L,n1)';
qw0 = qwf*ones(n1,1);
h0  = zeros(n1,1);
Fr1 = zeros(n1,1);
u1  = zeros(n1,1);
sf  = zeros(n1,1);
dh  = zeros(n1,1);

%Store: information of substrate stratigraphy;
%indup: number of sublayer at every node;
%Pup:proportion of sediment on uppermost sublayer;
%Lup: thickness of uppermost sublayer;
Store=cell(n,1);
indup=ceil((zb-La+4*Ls)./Ls);
indup(n_dam) = 1;
for is=1:n
    Store{is,1}=zeros(300,ng);
    Store{is,1}(1:indup(is),:)=ones(indup(is),1)*fbase;
end
fup=ones(n,1)*fbase;
Lup=zb+4*Ls-La-Ls*(indup-1);
Lup(n_dam) = Ls;
La(n_dam)  = 0;

%%%%%%Parameters in the Sediment Transport Module%%%%%%
%taub: bed shear stress;
%taur:reference shear stress;
%taurg:reference shear stress for surface geometric mean grain size;
%shirg: reference shields number for surface geometric mean grain size;
%ustar: shear velocity;
%b:component to compute taur;
%phi:taub on taur;
%Wstar: dimensionless bedload transport rate;
%qb: sediment transport rate;
%qbT: total volume sediment transport rate;
%qsT: total transport rate of sand;
%qgT: total transport rate of gravel;
%pbi: grain size distribution of bedload;
taub=zeros(n,1);
taur=zeros(n,ng);
taurg=zeros(n,1);
shirg=zeros(n,1);
ustar=zeros(n,1);
b=zeros(n,ng);
phi=zeros(n,ng);
Wstar=zeros(n,ng);
qb=zeros(n,ng);
qbT=zeros(n,1);
qsT=zeros(n,1);
qgT=zeros(n,1);
pbi=zeros(n,ng);
%fI: interfacial exchange fractions;
%dzb: change of bed evolution;
%dLa: change of activer layer thickness;
%delta: change of substrate elevation;
fI=zeros(n,ng);
dzb=zeros(n,1);
dLa=zeros(n,1);
dFi=zeros(n,ng);
delta=zeros(n,1);
%qwt: flow discharge per unit width at different time;
%zbt: bed elevation at different time;
%st: bed slope at different time;
%sct: (central scheme) bed slope at different time;
%Dsgt: Dsg at different time;
%FSt: FS at different time;
%Fit: Fi at different time;
%qbTt: total sediment transport rate per unit width at different time;
%Dg_loadt: Dg_load at different time;
qwt=zeros(n,Nt);
zbt=zeros(n,Nt);
dzbt=zeros(n,Nt);
st=zeros(n,Nt);
sct=zeros(n,Nt);
Dsgt=zeros(n,Nt);
qgTt=zeros(n,Nt);
qsTt=zeros(n,Nt);
taubt=zeros(n,Nt);
FRt=zeros(n,Nt);
FSt=zeros(n,Nt);
Fit=cell(Nt,1);
qbTt=zeros(n,Nt);
Dg_loadt=zeros(n,Nt);
it=1;


%%%%%%Temporal Evolution%%%%%%
while t<Duration
    %%%%%%Hydraulics£ºBackwater Equation%%%%%%
    zb1 = interp1(x, zb, x1);
    s1  = [(zb1(1:n1-1)-zb1(2:n1))/dx1; (zb1(n1-1)-zb1(n1))/dx1];
    D90s = interp1(x, D90, x1);
    hc=(qw0./0.75./g.^0.5).^(2/3);
    hn=(s1>0).*(qw0).^(3/5).*(nk.*D90s).^(1/10)./(alpr)^(3/5)./(g*s1).^(3/10)+(s1<=0).*hc;
    
    h0(n1)=hn(n1);
    for ih=1:n1-1
        Fr1(n1-ih+1)=qw0(n1-ih+1)/h0(n1-ih+1)./(g.*h0(n1-ih+1)).^0.5;
        if Fr1(n1-ih+1)<=0.75
            u1(n1-ih+1)=qw0(n1-ih+1)/h0(n1-ih+1);
            ustar(n1-ih+1)=u1(n1-ih+1)./alpr.*(nk.*D90s(n1-ih+1)./h0(n1-ih+1)).^(1/6);
            sf(n1-ih+1)=ustar(n1-ih+1).^2./g./h0(n1-ih+1);
            dh(n1-ih+1)=(s1(n1-ih)-sf(n1-ih+1))./(1-Fr1(n1-ih+1).^2);
            h0(n1-ih)=h0(n1-ih+1)-dh(n1-ih+1)*dx1;
        else
            h0(n1-ih)=hn(n1-ih);
        end
        Fr1(n1-ih)=qw0(n1-ih)/h0(n1-ih)./(g.*h0(n1-ih)).^0.5;
        if ((Fr1(n1-ih)<=0.75)&&(Fr1(n1-ih+1)>0.75))
            ha=hc(n1-ih+1);
            ua=qw0(n1-ih+1)/ha;
            Fra=qw0(n1-ih+1)/ha./(g.*ha).^0.5;
            ustara=ua./alpr.*(nk.*D90s(n1-ih+1)./ha).^(1/6);
            sfa=ustara.^2./g./ha;
            dha=(s1(n1-ih)-sfa)./(1-Fra.^2);
            h0(n1-ih)=ha-dha*dx1;
            Fr1(n1-ih)=qw0(n1-ih)/h0(n1-ih)./(g.*h0(n1-ih)).^0.5;
        end
        h0(n1-ih)=(Fr1(n1-ih)<=0.75)*h0(n1-ih)+(Fr1(n1-ih)>0.75)*hn(n1-ih);
        h0(n1-ih)=(h0(n1-ih)>=0.1)*h0(n1-ih)+(h0(n1-ih)<0.1)*0.1;
    end
    h = interp1(x1, h0, x);
    
    %u: flow velocity;
    %FR: Froude number;
    u=qw./h;
    FR=u./(g*h).^0.5;
    
    %%%%%%Sediment transport: Wilcock and Crowe (2003)%%%%%%
    Cf=(nk.*D90./h).^(1/3)/alpr^2;
    taub=rho.*Cf.*u.^2;
    ustar=(taub./rho).^0.5;
    shirg=0.021+0.015*exp(-20*FS);
    taurg=shirg.*R.*rho.*g.*Dsg;
    b=0.67./(1+exp(1.5-(ones(n,1)*D)./(Dsg*ones(1,ng))));
    taur=(taurg*ones(1,ng)).*((ones(n,1)*D)./(Dsg*ones(1,ng))).^b;
    phi=(taub*ones(1,ng))./taur;
    Wstar=(phi<1.35).*0.002.*(phi).^7.5+(phi>=1.35).*14.*(1-0.894./(phi).^0.5).^4.5;
    qb=Wstar.*Fi.*((ustar.^3)*ones(1,ng))./R./g;
    if etaa < La0
        qb(n_dam,:) = min([(1-p)*La(n_dam)*Fi(n_dam,:)*dx/dt/timeyear+qb(n_dam-1,:); qb(n_dam,:)]);
    end
    qbT=sum(qb,2);
    qsT=sum(qb(:,1:3),2);
    qgT=sum(qb(:,4:9),2);
    pbi=qb./(qbT*ones(1,ng));
    if qbT(n_dam) == 0
        pbi(n_dam,:) = 0;
    end
    
    %%%%%%Exner Equation%%%%%%
    %%%Evolution of Bed Elevation
    qbTback=[qtfT;qbT(1:n-1)];
    qbTit=qbT;
    qbTfrnt=[qbT(2:n);2*qbT(n)-qbT(n-1)];
    qbTdif=au*(qbTit-qbTback)+(1-au)*(qbTfrnt-qbTit);
    dzb=-qbTdif/dx/(1-p)*dt*timeyear;
    dzb(n)=0;
    zb=zb+dzb;
    %%%Evolution of Surface Fraction
    La_tmp = min([La0; etaa+dzb(n_dam)]);
    dLa(n_dam) = La_tmp-La(n_dam);
    delta=dzb-dLa;
    fI=((delta<=0)*ones(1,ng)).*fup+((delta>0)*ones(1,ng)).*(alpha*Fi+(1-alpha)*pbi);
    qbback=[qtf;qb(1:n-1,:)];
    qbit=qb;
    qbfrnt=[qb(2:n,:);2*qb(n,:)-qb(n-1,:)];
    qbdif=au*(qbit-qbback)+(1-au)*(qbfrnt-qbit);
    dFi=((-qbdif+qbTdif*ones(1,ng).*fI)/dx*dt*timeyear/(1-p)-(Fi-fI).*(dLa*ones(1,ng)))./(La*ones(1,ng));
    if etaa<La0 && etaa>0.05
        dFi(n_dam,:) = (-qbdif(n_dam,:)/dx*dt*timeyear/(1-p)+Fi(n_dam,:)*La(n_dam))/La_tmp-Fi(n_dam,:);
    elseif etaa<=0.05
        dFi(n_dam,:) = Fi(n_dam-1,:)-Fi(n_dam,:);
    end
    Fi=Fi+dFi;
    Fi=(Fi>0).*Fi;
    Fi=Fi./(sum(Fi,2)*ones(1,ng));
    La(n_dam) = La_tmp;
    etaa = etaa+dzb(n_dam);
    
    %%%%%%Stratigraphy Storage%%%%%%
    indn=find((delta<=Ls-Lup)&(delta>-Lup));
    indinc=find(delta>(Ls-Lup+ep));
    inddec=find(delta<=-Lup);
    fup(indn,:)=(fup(indn,:).*(Lup(indn)*ones(1,ng))+fI(indn,:).*(delta(indn)*ones(1,ng)))./((Lup(indn)+delta(indn))*ones(1,ng));
    Lup(indn)=Lup(indn)+delta(indn);
    if size(indinc,1)>0
        for is=1:size(indinc,1)
            ii=indinc(is);
            fup(ii,:)=(fup(ii,:).*Lup(ii)+fI(ii,:).*(Ls-Lup(ii)))./Ls;
            Store{ii,1}(indup(ii),:)=fup(ii,:);
            inc=ceil((delta(ii)-(Ls-Lup(ii)))/Ls);
            Store{ii,1}(indup(ii)+1:indup(ii)+inc,:)=ones(inc,1)*fI(ii,:);
            indup(ii)=indup(ii)+inc;
            fup(ii,:)=fI(ii,:);
            Lup(ii)=delta(ii)-(Ls-Lup(ii))-(inc-1)*Ls;
        end
    end
    if size(inddec,1)>0
        for is=1:size(inddec,1)
            id=inddec(is);
            dec=ceil((-Lup(id)-delta(id))/Ls);
            Store{id,1}(indup(id)-dec+1:indup(id),:)=zeros(dec,ng);
            indup(id)=indup(id)-dec;
            fup(id,:)=Store{id,1}(indup(id),:);
            Lup(id)=dec*Ls+Lup(id)+delta(id);
        end
    end
    
    %%%%%%Parameter Update%%%%%%
    t=t+dt;
    i=i+1;
    s=[(zb(1:n-1)-zb(2:n))/dx;(zb(n-1)-zb(n))/dx];
    FS=sum(Fi(:,1:3),2);
    Dsg=2.^(Fi*Phi')/1000;
    Cdg=[zeros(n,1),cumsum(Fi,2)];
    for iD=1:n
        ind=find(Cdg(iD,:)>=0.9,1);
        D1=Phin(ind-1);
        D2=Phin(ind);
        P1=Cdg(iD,ind-1);
        P2=Cdg(iD,ind);
        D90(iD,1)=2^((0.9-P1)/(P2-P1)*(D2-D1)+D1)/1000;
    end
    
    %%%Record Information Every Several Years
    if (mod(i,interval)==0)
        zbt(:,it)=zb;
        dzbt(:,it)=zb-zbini;
        st(:,it)=s;
        Dsgt(:,it)=Dsg;
        FSt(:,it)=FS;
        Fit{it,1}=Fi;
        qbTt(:,it)=qbT;
        qgTt(:,it)=qgT;
        qsTt(:,it)=qsT;
        it;
        it=it+1;
    end
end

for is=1:n
    Store{is,1}(indup(is),:)=fup(is,:);
end
toc