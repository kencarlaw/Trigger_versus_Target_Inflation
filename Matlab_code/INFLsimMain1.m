% INFLsimMain1.m 
% for "Monetary Policy in a Dynamic Model of Inflation"
% Kenneth I. Carlaw 21 Dec. 2024

clear

T=10000;
N=100;
NN=1;
TT=6;
gam=0.8;

%NN=N+1;
cr=0.012;
Istar=0.025;ci=0.15;b=1;%+cr;%+ci*Istar;
rho=0.5;
mu=1+0.000;sig=0.1;c=1;
q=1;qx=1.2;qxx=1;
adj=1;

p1=zeros(T,N,1);p1l=zeros(T,N,1);v=zeros(T,1);mp=zeros(T,1);
in=ones(T,1);in2=ones(T,1);win=ones(T,1);
ain=ones(T,1);
rh=zeros(T,1);fb=zeros(T,1);xin=zeros(T,1);
pr1=zeros(T,N,1);pr2=zeros(T,N,1);
hOm=zeros(T,N,1);hOml=zeros(T,N,1);hOme=zeros(T,N,1);hOmel=zeros(T,N,1);
hew=zeros(T,N,1);hewl=zeros(T,N,1);
hw=zeros(T,N,1);hwl=zeros(T,N,1);
mp1=zeros(T,1);mp1l=zeros(T,1);mhew=zeros(T,1);mhewl=zeros(T,1);test=zeros(T,1);
mhw=zeros(T,1);UE=zeros(T,1);


for t=2:T
    w=ones(TT,N,1);ew=ones(TT,N,1);
    p=ones(TT,N,1);
    infl=ones(TT,1);
    Om=zeros(TT,N,1);Ome=zeros(TT,N,1);
    if t>2
        for i=1:2
            p(i,:)=p1(t-1,:);
            ew(i,:)=hew(t-1,:);
            Ome(i,:)=hOme(t-1,:)*qxx;
        end
    else
        for i=1:5
            Ome(i,:)=adj*N-1;
            Om(i,:)=adj*N-1;
        end
    end
    ap=normrnd(mu,sig,TT,N,1);
    M=zeros(TT,1);EM=zeros(TT,1);
    epn=zeros(TT,N);exn=zeros(TT,N);pren=zeros(TT,N);
    ein=ones(TT,N);

    if in2(t-1)<0
        xin(t)=0.00000000001;
    else
        xin(t)=in2(t-1);
    end

    fb(t)=gam*min(1,(cr+ci*Istar)/(xin(t)));
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
        infl(m)=1+sum(log(p(m-1,:))-log(p(m-2,:)))/N;%sum(p(m-2,:)); 
        for k=1:N
            pb(k,:)=p(m-1,:);%(p(m-1,:)+p(m-2,:))/2;
            pa(k,:)=p(m-2,:);%(p(m-2,:)+p(m-3))/2;
        end
        for k=1:N
            for i=1:refsize
                for j=1:N
                    if j==ref(k,i)
                        pe(k,j)=pb(k,ref(k,i));
                        peh(k,j)=pa(k,j);
                        dpe(k,j)=log(pe(k,j))-log(pa(k,j));                    
                    end
                    %dpe(k,j)=pe(k,j)-pa(k,j);                    
                end
            end
            ein(m,k)=1+(sum(dpe(k,:))/(N-1));%*((N-eps)/eps);
            for j=1:N
                if pe(k,j)==1
                    pe(k,j)=pa(k,j);
                end
            end
        end

        for k=1:N
            for j=1:N-1
                if j~=k
                    Ome(m,k)=Ome(m,k)+(pe(k,j))^(rho/(rho-1));
                end
            end
            for j=1:N-1
                if j~=k
                    Om(m,k)=Om(m,k)+pb(k,j)^(rho/(rho-1));
                end
            end
        end
        EM(m)=sum(ew(m-1,:));
        for k=1:N
            ew(m,k)=(c/ap(m,k))*ew(m-1,k)*(ein(m-1,k)-qx*(infl(m-2)-ein(m-2,k)));
            w(m,k)=(c/ap(m,k))*ew(m-1,k)*(infl(m-1));
            %epn(m,k)=(ew(m,k))/(ap(m,k))*(1+sqrt(1+(ap(m,k))/(ew(m,k))*(1/(Ome(m,k)))));
            p(m,k)=(ew(m,k))/(1)*(1+sqrt(1+(1)/(ew(m,k))*(1/(Ome(m,k)))));
        end
        
    end
    in(t)=infl(TT-1)-1;
    p1(t,:)=p(TT,:);
    p1l(t,:)=p(TT-1,:);
    mp(t)=mean(p(TT,:));
    mp1l(t)=mean(p1l(t,:));
    hOm(t,:)=Om(TT,:);
    hOml(t,:)=Om(TT-1,:);
    hOme(t,:)=Ome(TT,:);
    hOmel(t,:)=Ome(TT-1,:);
    hewl(t,:)=ew(TT,:);
    hew(t,:)=ew(TT-1,:);
    hwl(t,:)=w(TT-1,:);
    hw(t,:)=w(TT,:);
    mhew(t)=mean(hew(t,:));
    mhewl(t)=mean(hewl(t,:));
    mhw(t)=mean(hwl(t,:));
    UE(t)=(1-(mhw(t)/mhew(t)));
    %test(t)=mean((ein(TT-1,:)-qx*(infl(TT-2)-ein(TT-2,:))));

    in(t)=sum(p1(t,:)-p1l(t,:))/sum(p1l(t,:));
    in2(t)=sum(log(p1(t,:))-log(p1l(t,:)))/N;
    win(t)=sum(abs(log(hew(t,:))-log(hwl(t,:))))/N;
    ain(t)=(in(t)+in(t-1))/2;
    for i=1:N
        if p1(t,i)>p1l(t,i)
            v(t)=v(t)+1;
        end
    end
    t
end

figure
histogram(in2(15:T))

figure
plot(rh(15:T))

figure
plot(in2(15:T))

figure
plot(UE(15:T))
hold on
plot(in2(15:T))
xlim([1 1000])
hold off

%save C:\Users\kcarlaw\Documents\MATLAB\INFL2024\p1 p1 -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\INFL2024\p1lag p1lag -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\INFL2024\w2 w2 -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\INFL2024\Omm Omm -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\INFL2024\ddp ddp -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\INFL2024\in in -ASCII -DOUBLE;

