%INFLmainconv.m
%Baseline "Monetary Policy in a Dynammic Model of Inflation"
%K. I. Carlaw Oct, 2024

clear

%parameters 

N=100;      %population of firms
NN=35;
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

%IND=[0.035 0.035];
IND=[0.0075 0.0575]; % Optimal two bin montery policy
Inst=0.0400;

critconv1=0.02;
%I=zeros(NN,1);
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
load iin2b26.mat
iin=iin2b26;
for k=1:1
    c=1;
    bb=0;
    %iin=zeros(NNN,NNN,1);
    Q=0;
    Pii=zeros(NNN,NNN,1);
    TMB=zeros(2,2);per=zeros(2,2);TMBC=zeros(2,2);
    MC=1;
    EDgb=zeros(MC,1);EDbb=zeros(MC,1);
    GBper=zeros(MC,1);BBper=zeros(MC,1);
    TGB=zeros(MC,1);TBB=zeros(MC,1);
    mciGA=zeros(MC,1);mciBA=zeros(MC,1);
    

    while Q < 1
        NCw=zeros(Block,1);I=zeros(Block,1);
        p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
        in2=ones(Block,1);
        ain=ones(Block,1);
        rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
        hOme=zeros(Block,N,1);
        hew=zeros(Block,N,1);    
        I(1)=IND(1);

        ip=zeros(Block,1);
        iGA=NaN(Block,1);iBA=NaN(Block,1);

        bb=bb+1
        
        for t=1:Block
            ew=ones(TT,N,1);
            p=ones(TT,N,1);
            infl=ones(TT,1);
            Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
            if t>2
                if in2(t-1) <= Inst(1)
                    I(t)=IND(1);
                else
                    I(t)=IND(2);
                end
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
        
            fb(t)=gam*min(1,(cr+ci*I(t))/(xin(t)));
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
            NCw(t)=lam2*I(t)+lam*in2(t);
            for i=1:NNN
                for j=1:NNN
                    if t>1
                        if (edges(i) <= in2(t-1)) && (in2(t-1) <= edges(i+1))
                            if (edges(j) <= in2(t)) && (in2(t) <= edges(j+1))
                                iin(i,j)=iin(i,j)+1;
                            end
                        end  
                    end
                end
            end
            if in2(t)<Inst(1)
                ip(t)=0;
                iGA(t)=in2(t);
            else
                ip(t)=1;
                iBA(t)=in2(t);
            end
            if t > 1
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
        if bb<2
            v=vw;
            NC=NCw;
            in=in2;
            rf=rh;
            iLA=iGA;
            iHA=iBA;
        else
            vhold=cat(2,v',vw');
            NChold=cat(2,NC',NCw');
            inhold=cat(2,in',in2');
            rfhold=cat(2,rf',rh');
            iLAhold=cat(2,iLA',iGA');
            iHAhold=cat(2,iHA',iBA');
            v=vhold';
            NC=NChold';
            in=inhold';
            rf=rfhold';
            iLA=iLAhold';
            iHA=iHAhold';
        end
        if bb>9999
            Q=1;
        end
        for j=1:2
            for i=1:2
                TMB(j,i)=TMBC(j,i)/sum(TMBC(j,:));
                per(j,i)=1/(1-TMB(j,i));
            end
        end
        EDgb(k)=per(1,1);
        EDbb(k)=per(2,2);
        GBper(k)=TMB(1,1);
        if isnan(GBper(k))
            GBper(k)=0;
        end
        BBper(k)=TMB(2,2);
        if isnan(BBper(k))
            BBper(k)=0;
        end    
        TGB(k)=(TMBC(1,1)+TMBC(1,2))/(bb*(Block-2));
        TBB(k)=(TMBC(2,2)+TMBC(2,1))/(bb*(Block-2));
        iLA(isnan(iLA))=[];
        iHA(isnan(iHA))=[];
        mciGA(k)=mean(iLA);
        mciBA(k)=mean(iHA);
        
    end
    for i=1:NNN
        for j=1:NNN
            Pii(i,j)=(iin(i,j)*10)/sum(iin(i,:)*10);
            if isnan(Pii(i,j))
                Pii(i,j)=0;
            end
        end
    end
    Eii=Pii*inp;
    U=Eii-inp;
    for i=1:NNN
        if U(i)>=0
            L(i)=1;
        else
            L(i)=-1;
        end
        if i>1            
            if L(i)-L(i-1)~=0
                count(i)=i-1;
            end
        end
    end
    count(count==0)=[];
    dEii=zeros(length(count),1);
    for i=1:length(count)
        dEii(i)=Eii(count(i));
    end
end
DEi=Eii-inp;
DIi=inp-inp;

xvec=zeros(NNN,2);
yvec=zeros(NNN,2);
for i=1:count(1)
    for j=1:2
        xvec(i,j)=0.033+(i-1)*0.005 + (j-1)*0.005;
        yvec(i,j)=0.621;
    end
end
for i=count(1)+1:NNN
    for j=1:2
        xvec(i,j)=0.089+NNN*0.009 -((i-1-count(1))*0.01 + (j-1)*0.01);
        yvec(i,j)=0.621;
    end
end
iinp=inp*100;
EEii=Eii*100;

figure
tile=tiledlayout(3,1);
tile.Padding='tight';
tile.TileSpacing='tight';
nexttile ([3 1])
hold on
box on
plot(inp*100,DEi*100,'Color','k','LineStyle','--','LineWidth',3)
plot(inp*100,DIi*100,'Color','k')
for i=1:(NNN)/3
    annotation('textarrow',xvec(i*3,:),yvec(i*3,:))
end 
xlabel('Inflation in previous period (\pi^{t-1})','FontSize',17)
ylabel('E(\pi^{t}|\pi^{t-1}) - \pi^{t}')
txt=['\pi^{t-1}=',num2str(EEii(count(1)),'%2.2f'),'%'];
text(iinp(count(1))-0.3,-0.03,txt,'FontSize',17);
xlim([0 20])
ylim([-0.9 0.6])
xtickformat('percentage')
ytickformat('percentage')
hold off


