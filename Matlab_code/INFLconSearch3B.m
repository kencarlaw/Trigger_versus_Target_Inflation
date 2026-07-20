%INFLmainconv.m
%Baseline "Monetary Policy in a Dynammic Model of Inflation"
%K. I. Carlaw July, 2026

clear

%parameters 

N=100;      %population of agents
NN=29;
NNN=N+1;
MM=500;

Block=1500;
TT=6;

    %baseline parameterization
gam=0.8;    % max objective apprehension prob.
ic=0.0025;
cr=0.012;
Istar=0;ci=0.15;
rho=0.5;
mu=1+0.000;sig=0.1;
qx=1.2;
lam=1;
lam2=1;

X=5;       %Set X <= 5

critconv1=0.03;

mv=zeros(NN,1);
mi=zeros(NN,1);
cost=zeros(NN,1);

edges=zeros(NNN+1,1);
for j=1:NNN+1
    edges(j)=(j-1)*0.0035-0.02;
end
MC=20;   % number of initial random policies. Run the program X times MC = 1000. If MC=20, X=1000/20

CM=zeros(MC,1);
IIIgb=zeros(MC,1);
IIIbb=zeros(MC,1);
IIItb=zeros(MC,1);
BnE=zeros(MC,1);
BnE2=zeros(MC,1);

for mc=1:MC
    %Initial values
    BE=ic*round(unifrnd(1,28));
    BE2=ic*round(unifrnd(round(BE/ic),28));
    Igb=ic*2;%ic*round(unifrnd(1,18));
    Itb=ic*round(unifrnd(round(Igb/ic),NN));
    Ibb=ic*round(unifrnd(18,NN));
    
    crit=1;
    IbbM=.1;
  
    while crit>0.5
        costold=IbbM;
    
        % Search for Bin 1 edge inflation
        BinE=zeros(2*X+1,1);
        Bcost=zeros(2*X+1,1);
        for k=1:2*X+1
            BinE(k)=0.02+((round(BE/ic)-1)+(k-1))*ic;
            c=1;
            bb=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            iconv=[];iconvlag=[];
            while ((cc<1) && (bb<MM))
                NCw=zeros(Block,1);
                p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
                in2=ones(Block,1);
                ain=ones(Block,1);
                rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
                hOme=zeros(Block,N,1);
                hew=zeros(Block,N,1);
                I=zeros(Block,1);
                   
                bb=bb+1;
                for t=1:Block
                    ew=ones(TT,N,1);
                    p=ones(TT,N,1);
                    infl=ones(TT,1);
                    Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
                    if bb>1
                        in2(1)=in(Block*(bb-1)-1);
                        in2(2)=in(Block*(bb-1)-1);
                        in2(3)=in(Block*(bb-1)-1);
                    end
                    I(1)=Igb;
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
                        if in2(t-1)<BinE(k)
                            I(t)=Igb;
                        elseif (BinE(k)<=in2(t-1)) && (in2(t-1)<BE2)
                            I(t)=Itb;
                        elseif BE2<= in2(t-1)
                            I(t)=Ibb;
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
                    % Best Response Price setting kernel
                    for m=3:TT
                        infl(m)=real(1+sum(log(p(m-1,:))-log(p(m-2,:)))/N); 
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
                            p(m,g)=real((ew(m,g))/(1)*(1+sqrt(1+(1)/(ew(m,g))*(1/(Ome(m,g))))));
                        end    
                    end
        
                    p1(t,:)=p(TT,:);
                    p1l(t,:)=p(TT-1,:);
                    mp(t)=mean(p(TT,:));
                    hOme(t,:)=Ome(TT,:);
                    hew(t,:)=ew(TT-1,:);
                    in2(t)=real(sum(log(p1(t,:))-log(p1l(t,:)))/N);
                    NCw(t)=lam2*I(t)+lam*in2(t);
                end
                if bb<2
                    NC=NCw;
                    in=in2;
                else
                    NChold=cat(2,NC',NCw');
                    inhold=cat(2,in',in2');
                    NC=NChold';
                    in=inhold';
                end
        
                iconv=in(1:bb*Block);
                iconvlag=in(1:(bb-1)*Block);
                freqiconv=histcounts(iconv(:),edges)/((bb)*Block);
                freqiconvlag=histcounts(iconvlag(:),edges)/((bb-1)*Block);
                TESTCONV1=zeros(1,NNN);
                if bb>1
                    for j=1:NNN
                        TESTCONV1(j)=abs(freqiconv(j)-freqiconvlag(j));
                    end
                    TESTCONV(bb)=sum(TESTCONV1);                
                end
                if bb > 4
                    if (TESTCONV(bb)<=critconv1) && (TESTCONV(bb-1)<=critconv1) ...
                        && (TESTCONV(bb-2)<=critconv1) && (TESTCONV(bb-3)<=critconv1)...
                        && (TESTCONV(bb-4)<=critconv1)
                        c=bb;
                        cc=1;
                    end
                end
            end
            Bcost(k)=mean(NC);  
            k
        end
        [BM,Ib]=min(Bcost);
        BE=BinE(Ib);





        % Search for Bin 2 edge inflation
        BinE2=zeros(2*X+1,1);
        Bcost2=zeros(2*X+1,1);
        for k=1:2*X+1
            if BE2<BE
                BE2=BE;
            end
            BinE2(k)=0.02+((round(BE2/ic)-1)+(k-1))*ic;
            c=1;
            bb=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            iconv=[];iconvlag=[];
            while ((cc<1) && (bb<MM))
                NCw=zeros(Block,1);
                p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
                in2=ones(Block,1);
                ain=ones(Block,1);
                rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
                hOme=zeros(Block,N,1);
                hew=zeros(Block,N,1);
                I=zeros(Block,1);
                   
                bb=bb+1;
                for t=1:Block
                    ew=ones(TT,N,1);
                    p=ones(TT,N,1);
                    infl=ones(TT,1);
                    Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
                    if bb>1
                        in2(1)=in(Block*(bb-1)-1);
                        in2(2)=in(Block*(bb-1)-1);
                        in2(3)=in(Block*(bb-1)-1);
                    end
                    I(1)=Igb;
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
                        if in2(t-1)<BE
                            I(t)=Igb;
                        elseif (BE<=in2(t-1)) && (in2(t-1)<BinE2(k))
                            I(t)=Itb;
                        elseif BinE2(k)<= in2(t-1)
                            I(t)=Ibb;
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
                    % Best Response Price setting kernel
                    for m=3:TT
                        infl(m)=real(1+sum(log(p(m-1,:))-log(p(m-2,:)))/N); 
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
                            p(m,g)=real((ew(m,g))/(1)*(1+sqrt(1+(1)/(ew(m,g))*(1/(Ome(m,g))))));
                        end    
                    end
        
                    p1(t,:)=p(TT,:);
                    p1l(t,:)=p(TT-1,:);
                    mp(t)=mean(p(TT,:));
                    hOme(t,:)=Ome(TT,:);
                    hew(t,:)=ew(TT-1,:);
                    in2(t)=real(sum(log(p1(t,:))-log(p1l(t,:)))/N);
                    NCw(t)=lam2*I(t)+lam*in2(t);
                end
                if bb<2
                    NC=NCw;
                    in=in2;
                else
                    NChold=cat(2,NC',NCw');
                    inhold=cat(2,in',in2');
                    NC=NChold';
                    in=inhold';
                end
        
                iconv=in(1:bb*Block);
                iconvlag=in(1:(bb-1)*Block);
                freqiconv=histcounts(iconv(:),edges)/((bb)*Block);
                freqiconvlag=histcounts(iconvlag(:),edges)/((bb-1)*Block);
                TESTCONV1=zeros(1,NNN);
                if bb>1
                    for j=1:NNN
                        TESTCONV1(j)=abs(freqiconv(j)-freqiconvlag(j));
                    end
                    TESTCONV(bb)=sum(TESTCONV1);                
                end
                if bb > 4
                    if (TESTCONV(bb)<=critconv1) && (TESTCONV(bb-1)<=critconv1) ...
                        && (TESTCONV(bb-2)<=critconv1) && (TESTCONV(bb-3)<=critconv1)...
                        && (TESTCONV(bb-4)<=critconv1)
                        c=bb;
                        cc=1;
                    end
                end
            end
            Bcost2(k)=mean(NC);  
            k
        end
        [BM,Ib]=min(Bcost2);
        BE2=BinE2(Ib);


        
        % Search for Igb
        IIgb=zeros(NN,1);
        Igbcost=zeros(NN,1);
        if Igb-X*ic<1*ic
            Igbd=1*ic;
        elseif Igb-X*ic>=NN*ic
            Igbd=(NN-1)*ic;
        else
            Igbd=Igb-X*ic;
        end
        if Igb+X*ic>=NN*ic
            Igbu=NN*ic;
        else
            Igbu=Igb+X*ic;
        end        
        for k=round(Igbd/ic):round(Igbu/ic)
            IIgb(k)=k*ic;    
            c=1;
            bb=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            iconv=[];iconvlag=[];
            while ((cc<1) && (bb<MM))
                NCw=zeros(Block,1);
                p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
                in2=ones(Block,1);
                ain=ones(Block,1);
                rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
                hOme=zeros(Block,N,1);
                hew=zeros(Block,N,1);
                I=zeros(Block,1);
                   
                bb=bb+1;
                for t=1:Block
                    ew=ones(TT,N,1);
                    p=ones(TT,N,1);
                    infl=ones(TT,1);
                    Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
                    if bb>1
                        in2(1)=in(Block*(bb-1)-1);
                        in2(2)=in(Block*(bb-1)-1);
                        in2(3)=in(Block*(bb-1)-1);
                        %p(1,:)=ph(:);
                        %p(2,:)=ph(:);
                    end
                    I(1)=Igb;
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
                        if in2(t-1)<BE
                            I(t)=Igb;
                        elseif (BE<=in2(t-1)) && (in2(t-1)<BE2)
                            I(t)=Itb;
                        elseif BE2<= in2(t-1)
                            I(t)=Ibb;
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
                    % Best Response Price setting kernel
                    for m=3:TT
                        infl(m)=real(1+sum(log(p(m-1,:))-log(p(m-2,:)))/N); 
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
                            p(m,g)=real((ew(m,g))/(1)*(1+sqrt(1+(1)/(ew(m,g))*(1/(Ome(m,g))))));
                        end    
                    end
        
                    p1(t,:)=p(TT,:);
                    p1l(t,:)=p(TT-1,:);
                    mp(t)=mean(p(TT,:));
                    hOme(t,:)=Ome(TT,:);
                    hew(t,:)=ew(TT-1,:);
                    in2(t)=real(sum(log(p1(t,:))-log(p1l(t,:)))/N);
                    NCw(t)=lam2*I(t)+lam*in2(t);
                end
                if bb<2
                    NC=NCw;
                    in=in2;
                else
                    NChold=cat(2,NC',NCw');
                    inhold=cat(2,in',in2');
                    NC=NChold';
                    in=inhold';
                end
        
                iconv=in(1:bb*Block);
                iconvlag=in(1:(bb-1)*Block);
                freqiconv=histcounts(iconv(:),edges)/((bb)*Block);
                freqiconvlag=histcounts(iconvlag(:),edges)/((bb-1)*Block);
                TESTCONV1=zeros(1,NNN);
                if bb>1
                    for j=1:NNN
                        TESTCONV1(j)=abs(freqiconv(j)-freqiconvlag(j));
                    end
                    TESTCONV(bb)=sum(TESTCONV1);                
                end
                if bb > 4
                    if (TESTCONV(bb)<=critconv1) && (TESTCONV(bb-1)<=critconv1) ...
                        && (TESTCONV(bb-2)<=critconv1) && (TESTCONV(bb-3)<=critconv1)...
                        && (TESTCONV(bb-4)<=critconv1)
                        c=bb;
                        cc=1;
                    end
                end
                ph=p1(t,:);
            end
            Igbcost(k)=mean(NC);  
            k
        end
        [IgbM, Xgb]=min(Igbcost(round(Igbd/ic):round(Igbu/ic)));
        Igb=(Xgb-1)*ic+Igbd;
        


        % Search for Itb
        IItb=zeros(NN,1);
        Itbcost=zeros(NN,1);
        if Itb-X*ic<Igb
            Itbd=Igb;
        elseif Itb-X*ic>=NN*ic
            Itbd=(NN-1)*ic;
        else
            Itbd=Itb-X*ic;
        end
        if Itb+X*ic>=NN*ic
            Itbu=(NN-1)*ic;
        else
            Itbu=Itb+X*ic;
        end        
        if Itbu>Ibb
            Itbu=Ibb;
        end
        if Itbd > Itbu
            Itbd=Itbu;
        end
        for k=round(Itbd/ic):round(Itbu/ic)
            IItb(k)=k*ic;    
            c=1;
            bb=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            iconv=[];iconvlag=[];
            while ((cc<1) && (bb<MM))
                NCw=zeros(Block,1);
                p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
                in2=ones(Block,1);
                ain=ones(Block,1);
                rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
                hOme=zeros(Block,N,1);
                hew=zeros(Block,N,1);
                I=zeros(Block,1);
                   
                bb=bb+1;
                for t=1:Block
                    ew=ones(TT,N,1);
                    p=ones(TT,N,1);
                    infl=ones(TT,1);
                    Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
                    if bb>1
                        in2(1)=in(Block*(bb-1)-1);
                        in2(2)=in(Block*(bb-1)-1);
                        in2(3)=in(Block*(bb-1)-1);
                        %p(1,:)=ph(:);
                        %p(2,:)=ph(:);
                    end
                    I(1)=Igb;
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
                        if in2(t-1)<BE
                            I(t)=Igb;
                        elseif (BE<=in2(t-1)) && (in2(t-1)<BE2)
                            I(t)=Itb;
                        elseif BE2<= in2(t-1)
                            I(t)=Ibb;
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
                    % Best Response Price setting kernel
                    for m=3:TT
                        infl(m)=real(1+sum(log(p(m-1,:))-log(p(m-2,:)))/N); 
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
                            p(m,g)=real((ew(m,g))/(1)*(1+sqrt(1+(1)/(ew(m,g))*(1/(Ome(m,g))))));
                        end    
                    end
        
                    p1(t,:)=p(TT,:);
                    p1l(t,:)=p(TT-1,:);
                    mp(t)=mean(p(TT,:));
                    hOme(t,:)=Ome(TT,:);
                    hew(t,:)=ew(TT-1,:);
                    in2(t)=real(sum(log(p1(t,:))-log(p1l(t,:)))/N);
                    NCw(t)=lam2*I(t)+lam*in2(t);
                end
                if bb<2
                    NC=NCw;
                    in=in2;
                else
                    NChold=cat(2,NC',NCw');
                    inhold=cat(2,in',in2');
                    NC=NChold';
                    in=inhold';
                end
        
                iconv=in(1:bb*Block);
                iconvlag=in(1:(bb-1)*Block);
                freqiconv=histcounts(iconv(:),edges)/((bb)*Block);
                freqiconvlag=histcounts(iconvlag(:),edges)/((bb-1)*Block);
                TESTCONV1=zeros(1,NNN);
                if bb>1
                    for j=1:NNN
                        TESTCONV1(j)=abs(freqiconv(j)-freqiconvlag(j));
                    end
                    TESTCONV(bb)=sum(TESTCONV1);                
                end
                if bb > 4
                    if (TESTCONV(bb)<=critconv1) && (TESTCONV(bb-1)<=critconv1) ...
                        && (TESTCONV(bb-2)<=critconv1) && (TESTCONV(bb-3)<=critconv1)...
                        && (TESTCONV(bb-4)<=critconv1)
                        c=bb;
                        cc=1;
                    end
                end
                ph=p1(t,:);
            end
            Itbcost(k)=mean(NC);  
            k
        end
        [ItbM, Xtb]=min(Itbcost(round(Itbd/ic):round(Itbu/ic)));
        Itb=(Xtb-1)*ic+Itbd;




        
        % Search for Ibb
        IIbb=zeros(NN,1);
        Ibbcost=zeros(NN,1);
        if Ibb<Itb
            Ibb=Itb;
        end
        if Ibb-X*ic<1*ic
            Ibbd=1*ic;
        elseif Ibb-X*ic>=NN*ic
            Ibbd=(NN-1)*ic;
        else
            Ibbd=Ibb-X*ic;
        end
        if Ibb+X*ic>=NN*ic
            Ibbu=NN*ic;
        else
            Ibbu=Ibb+X*ic;
        end        
        for k=round(Ibbd/ic):round(Ibbu/ic)
            IIbb(k)=k*ic;    
            c=1;
            bb=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            iconv=[];iconvlag=[];
            while ((cc<1) && (bb<MM))
                NCw=zeros(Block,1);
                p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
                in2=ones(Block,1);
                ain=ones(Block,1);
                rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
                hOme=zeros(Block,N,1);
                hew=zeros(Block,N,1);
                I=zeros(Block,1);
                   
                bb=bb+1;
                for t=1:Block
                    ew=ones(TT,N,1);
                    p=ones(TT,N,1);
                    infl=ones(TT,1);
                    Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
                    if bb>1
                        in2(1)=in(Block*(bb-1)-1);
                        in2(2)=in(Block*(bb-1)-1);
                        in2(3)=in(Block*(bb-1)-1);
                        %p(1,:)=ph(:);
                        %p(2,:)=ph(:);
                    end
                    I(1)=Igb;
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
                        if in2(t-1)<BE
                            I(t)=Igb;
                        elseif (BE<=in2(t-1)) && (in2(t-1)<BE2)
                            I(t)=Itb;
                        elseif BE2<= in2(t-1)
                            I(t)=Ibb;
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
                    % Best Response Price setting kernel
                    for m=3:TT
                        infl(m)=real(1+sum(log(p(m-1,:))-log(p(m-2,:)))/N); 
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
                            p(m,g)=real((ew(m,g))/(1)*(1+sqrt(1+(1)/(ew(m,g))*(1/(Ome(m,g))))));
                        end    
                    end
        
                    p1(t,:)=p(TT,:);
                    p1l(t,:)=p(TT-1,:);
                    mp(t)=mean(p(TT,:));
                    hOme(t,:)=Ome(TT,:);
                    hew(t,:)=ew(TT-1,:);
                    in2(t)=real(sum(log(p1(t,:))-log(p1l(t,:)))/N);
                    NCw(t)=lam2*I(t)+lam*in2(t);
                end
                if bb<2
                    NC=NCw;
                    in=in2;
                else
                    NChold=cat(2,NC',NCw');
                    inhold=cat(2,in',in2');
                    NC=NChold';
                    in=inhold';
                end
        
                iconv=in(1:bb*Block);
                iconvlag=in(1:(bb-1)*Block);
                freqiconv=histcounts(iconv(:),edges)/((bb)*Block);
                freqiconvlag=histcounts(iconvlag(:),edges)/((bb-1)*Block);
                TESTCONV1=zeros(1,NNN);
                if bb>1
                    for j=1:NNN
                        TESTCONV1(j)=abs(freqiconv(j)-freqiconvlag(j));
                    end
                    TESTCONV(bb)=sum(TESTCONV1);                
                end
                if bb > 4
                    if (TESTCONV(bb)<=critconv1) && (TESTCONV(bb-1)<=critconv1) ...
                        && (TESTCONV(bb-2)<=critconv1) && (TESTCONV(bb-3)<=critconv1)...
                        && (TESTCONV(bb-4)<=critconv1)
                        c=bb;
                        cc=1;
                    end
                end
                ph=p1(t,:);
            end
            Ibbcost(k)=mean(NC);  
            k
        end
        [IbbM, Xbb]=min(Ibbcost(round(Ibbd/ic):round(Ibbu/ic)));
        Ibb=(Xbb-1)*ic+Ibbd;
        if Ibb<Igb
            Ibb=Igb;
        end        
        crit=abs(costold-IbbM)/IbbM;
    end
    BnE(mc)=BE;
    BnE2(mc)=BE2;
    IIIgb(mc)=Igb;
    IIItb(mc)=Itb;
    IIIbb(mc)=Ibb;
    CM(mc)=IbbM;
    
    mc
end

opt=[BE BE2 Igb Ibb Itb IbbM];

OPTB=[BnE BnE2 IIIgb IIItb IIIbb CM];
