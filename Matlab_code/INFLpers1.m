%INFLmainconv.m
%Baseline "Monetary Policy in a Dynammic Model of Inflation"
%K. I. Carlaw Oct, 2024

clear

%parameters 

N=100;      %population of firms
%NN=29;
NNN=N+1;
MM=500;

Block=10000;
TT=6;

    %baseline parameterization
gam=0.8;    % max infprmation prob. weight
ic=0.0025;
cr=0.012;
Istar=0;ci=0.15;
rho=0.5;
mu=1+0.000;sig=0.1;
qx=1.2;
lam=1;
lam2=1;

BE=[0.0865 0.0876 0.0886 0.0896 0.0924 0.0962 0.0988 0.1096 0.1154 0.1233 0.1244 0.1306 0.1356 0.1375 0.1397 0.1410 0.1419];
IND=[5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21];   % montery policy index 
NN=length(IND);

critconv1=0.02;
I=zeros(NN,1);
mv=zeros(NN,1);
mi=zeros(NN,1);
cost=zeros(NN,1);
icon=zeros(NN,5*Block,1);rfcon=zeros(NN,5*Block,1);

edges=zeros(NNN+1,1);
for j=1:NNN+1
    edges(j)=(j-1)*0.0035-0.02;
end

Eii=zeros(NNN,1);
U=zeros(NNN,1);
L=zeros(NNN,1);
count=zeros(NNN,1);
inp=zeros(NNN,1);
for i=1:NNN
    inp(i)=edges(i+1);
end

MC=25;

mGBper=zeros(NN,1);
mBBper=zeros(NN,1);
mTGB=zeros(NN,1);mTGB2=zeros(NN,1);
mTBB=zeros(NN,1);mTBB2=zeros(NN,1);
miGA=zeros(NN,1);
miBA=zeros(NN,1);

for k=1:NN
    I(k)=Istar+(IND(k)-1)*ic;
    c=1;
    iin=zeros(NNN,NNN,1);
    Pii=zeros(NNN,NNN,1);
    EDgb=zeros(MC,1);EDbb=zeros(MC,1);
    GBper=zeros(MC,1);BBper=zeros(MC,1);
    TGB=zeros(MC,1);TBB=zeros(MC,1);
    mciGA=zeros(MC,1);mciBA=zeros(MC,1);

    for mc=1:MC
        NCw=zeros(Block,1);
        p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
        in2=ones(Block,1);
        ain=ones(Block,1);
        rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
        hOme=zeros(Block,N,1);
        hew=zeros(Block,N,1);

        ip=zeros(Block,1);TMBC=zeros(2,2);
        iGA=NaN(Block,1);iBA=NaN(Block,1);
            
        for t=1:Block
            ew=ones(TT,N,1);
            p=ones(TT,N,1);
            infl=ones(TT,1);
            Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
            if t>2
                for i=1:2
                    p(i,:)=p1(t-1,:);
                    ew(i,:)=hew(t-1,:);
                    Ome(i,:)=hOme(t-1,:);
                end
                if in2(t-1)<0
                    xin(t)=0.00000000001;
                else
                    xin(t)=in2(t-1);
                end                
            else
                for i=1:5
                    Ome(i,:)=N-1;
                    Om(i,:)=N-1;
                end
            end
            ap=normrnd(mu,sig,TT,N,1);
            M=zeros(TT,1);EM=zeros(TT,1);
            epn=zeros(TT,N);exn=zeros(TT,N);pren=zeros(TT,N);
            pro1=zeros(TT,N);pro2=zeros(TT,N);
            ein=ones(TT,N);
        
            fb(t)=gam*min(1,(cr+ci*I(k))/(xin(t)));
            %bb(t)=cr+ci*I(k);
            refsize=binornd(N-1,fb(t));
            %refsize=round(N-(1-gam*min(1,fb(t)))*eps);
            %refsize=round(N-(1-gam*(1-1/(ep^(Istar/(in(t-1))))))*eps);
            %refsize=N-eps;
            if refsize < 0
                refsize=0;
            end
            rh(t)=refsize;
            
            refg=zeros(N,refsize+1);ref=zeros(N,refsize);
            Nvec=zeros(N,1);
            Nvec2=zeros(N,1);
            for i=1:N
                Nvec(i)=i;
            end 
            for i=1:N
                Nvec2=Nvec;
                Nvec2(i)=[];
                ref(i,:)=randsample(Nvec2,refsize);
            end
            pe=ones(N,N,1);pa=ones(N,N,1);pb=ones(N,N,1);
            peh=ones(N,N,1);dpe=zeros(N,N,1);

            for m=3:TT
                infl(m)=1+sum(log(p(m-1,:))-log(p(m-2,:)))/N; 
                for i=1:N
                    pb(i,:)=p(m-1,:);
                    pa(i,:)=p(m-2,:);
                end
                for g=1:N
                    for i=1:refsize
                        for j=1:N
                            if j==ref(g,i)
                                pe(g,j)=pb(g,ref(g,i));
                                peh(g,j)=pa(g,j);
                                dpe(g,j)=log(pe(g,j))-log(pa(g,j));                    
                            end
                        end
                    end
                    ein(m,g)=1+(sum(dpe(g,:))/(N-1));
                    for j=1:N
                        if pe(g,j)==1
                            pe(g,j)=pa(g,j);
                        end
                    end
                end
        
                for g=1:N
                    for j=1:N-1
                        if j~=g
                            Ome(m,g)=Ome(m,g)+(pe(g,j))^(rho/(rho-1));
                        end
                    end
                end
                EM(m)=sum(ew(m-1,:));
                for g=1:N
                    ew(m,g)=(c/ap(m,g))*ew(m-1,g)*(ein(m-1,g)-qx*(infl(m-2)-ein(m-2,g)));
                    p(m,g)=(ew(m,g))/(1)*(1+sqrt(1+(1)/(ew(m,g))*(1/(Ome(m,g)))));
                end    
            end

            p1(t,:)=p(TT,:);
            p1l(t,:)=p(TT-1,:);
            mp(t)=mean(p(TT,:));
            hOme(t,:)=Ome(TT,:);
            hew(t,:)=ew(TT-1,:);
            in2(t)=sum(log(p1(t,:))-log(p1l(t,:)))/N;
            for i=1:N
                if p1(t,i)>p1l(t,i)
                    vw(t)=vw(t)+1;
                end
            end        
            %NCw(t)=lam2*I(k)+lam*in2(t);
            if t>2  
                if in2(t)<BE(k)
                    ip(t)=0;
                    iGA(t)=in2(t);
                else
                    ip(t)=1;
                    iBA(t)=in2(t);
                end
                if (ip(t)==0) && (ip(t-1)==0)
                    TMBC(1,1)=TMBC(1,1)+1;
                elseif (ip(t)==0) && (ip(t-1)==1)
                    TMBC(1,2)=TMBC(1,2)+1;
                elseif (ip(t)==1) && (ip(t-1)==1)
                    TMBC(2,2)=TMBC(2,2)+1;
                else
                    TMBC(2,1)=TMBC(2,1)+1;
                end
            end
        end
        TMB=zeros(2,2);per=zeros(2,2);
        for j=1:2
            for i=1:2
                TMB(j,i)=TMBC(j,i)/sum(TMBC(j,:));
                per(j,i)=1/(1-TMB(j,i));
            end
        end
        EDgb(mc)=per(1,1);
        EDbb(mc)=per(2,2);
        GBper(mc)=TMB(1,1);
        if isnan(GBper(mc))
            GBper(mc)=0;
        end
        BBper(mc)=TMB(2,2);
        if isnan(BBper(mc))
            BBper(mc)=0;
        end    
        TGB(mc)=(TMBC(1,1)+TMBC(1,2))/(Block-2);
        TBB(mc)=(TMBC(2,2)+TMBC(2,1))/(Block-2);
        iGA(isnan(iGA))=[];
        iBA(isnan(iBA))=[];
        mciGA(mc)=mean(iGA);
        mciBA(mc)=mean(iBA);

        %if bb<2
        %    v=vw;
        %    %NC=NCw;
        %    in=in2;
        %    rf=rh;
        %else
        %    vhold=cat(2,v',vw');
        %    %NChold=cat(2,NC',NCw');
        %    inhold=cat(2,in',in2');
        %    rfhold=cat(2,rf',rh');
        %    v=vhold';
        %    %NC=NChold';
        %    in=inhold';
        %    rf=rfhold';
        %end
    end
    GBper(GBper==0)=[];
    BBper(BBper==0)=[];
    TGB(TGB==0)=[];
    TBB(TBB==0)=[];
    mciGA(mciGA==0)=[];
    mciGA(isnan(mciGA))=[];
    mciBA(mciBA==0)=[];
    mciBA(isnan(mciBA))=[];
    mGBper(k)=mean(GBper);
    BBper(isnan(BBper))=[];
    mBBper(k)=mean(BBper);
    mTGB(k)=mean(TGB);
    mTBB(k)=mean(TBB);
    miGA(k)=mean(mciGA);
    miBA(k)=mean(mciBA);     
    I(k)
end

load mibl26.mat
Ip=zeros(NN,1);
for i=1:NN
    Ip(i)=I(i)*100;
end
figure
tile=tiledlayout(4,5);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([1 5])
hold on
box on
%set(gca,'TickDir','out','TickLength',[0.005,0.005],'FontSize',11)
plot(Ip,mibl(5:20)*100,'Color','k','LineStyle','-','LineWidth',1.5)
plot(Ip,miGA*100,'Color','k','LineStyle',':','Marker','square','MarkerIndices',1:3:NN,'LineWidth',1.5)
plot(Ip,miBA*100,'Color','k','LineStyle','--','Marker','o','MarkerIndices',1:3:NN,'LineWidth',1.5)
xlim([0.01*100 0.0475*100])
ylim([0*100 0.23*100])
xtickformat('percentage')
ytickformat('percentage')
title('Panel 1: E(\pi|i) and FPs within the q-attractors','FontSize',20,'FontWeight','normal')
%xlabel('Enforcement resources (R)')
%ylabel('Expected violations','FontSize',14)
legend('Expected \pi','Focal Point in LA','Focal Point in HA','FontSize',15,'Location','west')
legend boxoff
legend('AutoUpdate','off')
%xline(55)
%yline(105)
hold off
nexttile([2 5])
hold on
%set(gca,'TickDir','out','TickLength',[0.005,0.005],'FontSize',11)
box on
plot(Ip,mGBper,'Color','k','LineStyle',':','Marker','square','MarkerIndices',1:3:NN,'LineWidth',1.5)
plot(Ip,mBBper,'Color','k','LineStyle','--','Marker','o','MarkerIndices',1:3:NN,'LineWidth',1.5)
title('Panel 2: Persistence of the q-attractors','FontSize',20,'FontWeight','normal')
legend('Persistence of LA','Persistence of HA','FontSize',15,'Location','west')
legend boxoff
%xlim([39 43])
%xticks([39 40 41 42 43])
ylim([0.78 1.01])
xlim([0.01*100 0.0475*100])
xtickformat('percentage')
%xticks([20 23 27 30 34 37 38 39 40 41 42 43 44])
legend('AutoUpdate','off')
%xline(55)
%yline(1.01)
hold off
nexttile([1 5])
hold on
%set(gca,'TickDir','out','TickLength',[0.005,0.005],'FontSize',11)
box on
plot(Ip,mTGB*100,'Color','k','LineStyle',':','Marker','square','MarkerIndices',1:3:NN,'LineWidth',1.5)
plot(Ip,mTBB*100,'Color','k','LineStyle','--','Marker','o','MarkerIndices',1:3:NN,'LineWidth',1.5)
title('Panel 3: % of time spent in q-attractors','FontSize',20,'FontWeight','normal') 
legend('% time in LA','% time in HA','FontSize',15,'Location','west')
legend boxoff
xlabel('Monetary policy (i)','FontSize',17)
%xlim([0 100])
%xticks([20 23 27 30 34 37 38 39 40 41 42 43 44])
ylim([-2 102])
xlim([0.01*100 0.0475*100])
xtickformat('percentage')
ytickformat('percentage')
legend('AutoUpdate','off')
%xline(55)
%yline(105)
hold off





