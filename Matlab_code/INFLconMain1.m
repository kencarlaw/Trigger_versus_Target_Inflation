%INFLmainconv.m
%Baseline "Monetary Policy in a Dynammic Model of Inflation"
%K. I. Carlaw July, 2026

clear

%parameters 

N=100;      %population of firms
NN=29;
NNN=N+1;
MM=500;

Block=20000;
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

critconv1=0.001;
I=zeros(NN,1);II=zeros(NN,1);
mv=zeros(NN,1);
mi=zeros(NN,1);
msi=zeros(NN,1);
cost=zeros(NN,1);
icon=zeros(NN,5*Block,1);rfcon=zeros(NN,5*Block,1);
qin=zeros(NN,1);qqin=zeros(NN,1);qn=zeros(NN,1);

edges=zeros(NNN+1,1);
for j=1:NNN+1
    edges(j)=(j-1)*0.0035-0.02;
end

for k=1:NN
    I(k)=Istar+(k-1)*ic;
    %II(k)=1+cr;%+I(k)*ci;
    c=1;
    bb=0;
    cc=0;
    TESTCONV=zeros(MM,1);
    TESTCONV(1)=1;
    while ((cc<1) && (bb<MM))
        NCw=zeros(Block,1);
        p1=ones(Block,N,1);p1l=ones(Block,N,1);vw=zeros(Block,1);mp=zeros(Block,1);
        in2=ones(Block,1);mein=zeros(Block,1);
        ain=ones(Block,1);
        rh=zeros(Block,1);fb=zeros(Block,1);xin=ones(Block,1);
        hOme=zeros(Block,N,1);
        hew=zeros(Block,N,1);    
            
        bb=bb+1;
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
            mein(t)=mean((ein(TT-1,:)-qx*(infl(TT-2)-ein(TT-2,:))))-1;
            for i=1:N
                if p1(t,i)>p1l(t,i)
                    vw(t)=vw(t)+1;
                end
            end        
            NCw(t)=lam2*I(k)+lam*in2(t);
        end
        if bb<2
            v=vw;
            NC=NCw;
            in=in2;
            rf=rh;
            mmein=mein;
        else
            vhold=cat(2,v',vw');
            NChold=cat(2,NC',NCw');
            inhold=cat(2,in',in2');
            rfhold=cat(2,rf',rh');
            mmeinhold=cat(2,mein',mmein');
            v=vhold';
            NC=NChold';
            in=inhold';
            rf=rfhold';
            mmein=mmeinhold';
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
    icon(k,:)=in(1:5*Block);
    rfcon(k,:)=rf(1:5*Block);
    mi(k)=mean(in);
    msi(k)=mean(mmein);
    mv(k)=mean(v);
    cost(k)=mean(NC);  
    qin(k)=sum(in>=-0.02);
    qqin(k)=sum(in<=.3335);
    qn(k)=length(in);
    I(k)
end


%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2023postRED\RRC.mat','RRC')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2023postRED\NRRC.mat','NRRC')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2023postRED\CMPL.mat','CMPL')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASBfinal2023\sdvbl.mat','sdvbl')

[acf38,lag38]=autocorr(icon(5,4:5*Block));
[acf39,lag39]=autocorr(icon(10,4:5*Block));
[acf40,lag40]=autocorr(icon(11,4:5*Block));
[acf41,lag41]=autocorr(icon(18,4:5*Block));
[acf42,lag42]=autocorr(icon(23,4:5*Block));
[acf43,lag43]=autocorr(icon(22,4:5*Block));

figure %figure 3.1
tile=tiledlayout(10,5);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile(1,[1 3])
plot(icon(5,1000:6000)*100,'Color',[0.5 0.5 0.5])
ylabel('{\it{i}} = 1.00%','FontSize',15)
ytickformat('percentage');
ylim([-0.02*100 0.3*100])
xlim([1 5000])
nexttile(6,[1 3]) 
plot(rfcon(5,1000:6000),'Linestyle','--','Color',[0.5 0.5 0.5])
ylabel('Info','FontSize',15)
ylim([0 100])
xlim([1 5000])
nexttile(4,[2 1])
histogram(icon(5,1000:6000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',0.01)
xlim([-0.02 0.3])
ylim([0 0.15])
view(90,-90)
nexttile(5,[2 1])
stem(lag38(2:4),acf38(2:4),'Filled','Color','k')
xlim([0 4]);
%title('Autocorrelation in violations')
%xlabel('Lags')
nexttile(11,[1 3])
plot(icon(10,1000:6000)*100,'Color',[0.5 0.5 0.5])
ylabel('{\it{i}} = 2.25%','FontSize',15)
ytickformat('percentage');
ylim([-0.02*100 0.3*100])
xlim([1 5000])
nexttile(16,[1 3]) 
plot(rfcon(10,1000:6000),'Linestyle','--','Color',[0.5 0.5 0.5])
ylabel('Info','FontSize',15)
ylim([0 100])
xlim([1 5000])
nexttile(14,[2 1])
histogram(icon(10,1000:6000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',0.01)
xlim([-0.02 0.3])
ylim([0 0.15])
view(90,-90)
nexttile(15,[2 1])
stem(lag39(2:4),acf39(2:4),'Filled','Color','k')
xlim([0 4]);
%title('Autocorrelation in violations')
%xlabel('Lags')
nexttile(21,[1 3])
plot(icon(11,1000:6000)*100,'Color',[0.5 0.5 0.5])
ylabel('{\it{i}} = 2.50%','FontSize',15)
ytickformat('percentage');
ylim([-0.02*100 0.3*100])
xlim([1 5000])
nexttile(26,[1 3]) 
plot(rfcon(11,1000:6000),'Linestyle','--','Color',[0.5 0.5 0.5])
ylabel('Info','FontSize',15)
ylim([0 100])
xlim([1 5000])
nexttile(24,[2 1])
histogram(icon(11,1000:6000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',0.01)
xlim([-0.02 0.3])
ylim([0 0.15])
view(90,-90)
nexttile(25,[2 1])
stem(lag40(2:4),acf40(2:4),'Filled','Color','k')
xlim([0 4]);
%title('Autocorrelation in violations')
%xlabel('Lags')
nexttile(31,[1 3])
plot(icon(18,1000:6000)*100,'Color',[0.5 0.5 0.5])
ylabel('{\it{i}} = 4.25%','FontSize',15)
ytickformat('percentage');
ylim([-0.02*100 0.3*100])
xlim([1 5000])
nexttile(36,[1 3]) 
plot(rfcon(18,1000:6000),'Linestyle','--','Color',[0.5 0.5 0.5])
ylabel('Info','FontSize',15)
ylim([0 100])
xlim([1 5000])
nexttile(34,[2 1])
histogram(icon(18,1000:6000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',0.01)
xlim([-0.02 0.3])
ylim([0 0.15])
view(90,-90)
nexttile(35,[2 1])
stem(lag41(2:4),acf41(2:4),'Filled','Color','k')
xlim([0 4]);
nexttile(41,[1 3])
plot(icon(23,1000:6000)*100,'Color',[0.5 0.5 0.5])
%title('Panel 6: R = 43')
ytickformat('percentage');
ylim([-0.01*100 0.3*100])
xlim([1 5000])
ylabel('{\it{i}} = 5.50%','FontSize',15)
nexttile(46,[1 3]) 
plot(rfcon(23,1000:6000),'Linestyle','--','Color',[0.5 0.5 0.5])
ylabel('Info','FontSize',15)
ylim([0 100])
xlim([1 5000])
xlabel('Period (T = 0,...,5000)','FontSize',20)
nexttile(44,[2 1])
histogram(icon(23,1000:6000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',0.01)
xlim([-0.01 0.3])
ylim([0 0.15])
ylabel('Inflation frequency','FontSize',20)
view(90,-90)
nexttile(45,[2 1])
stem(lag42(2:4),acf42(2:4),'Filled','Color','k')
xlim([0 4]);
ylim([0 1]);
%title('Autocorrelation in violations')
xlabel('AC for 3 lags','FontSize',20)


%figure % Figure 3.6
%tile=tiledlayout(3,1);
%tile.Padding='none';
%tile.TileSpacing='tight';
%nexttile ([3 1])
%hold on
%plot(I,mv,'Color','k','LineStyle','-','LineWidth',0.5)
%box on
%ylabel('Number of price changes')
%title('Panel 1: Expected price changes in the estimated stationary distribution')
%legend('Mean price changes E(\Delta p|I), D-model')
%legend boxoff
%xlabel('Interest rate (I)')
%hold off

load mibl26.mat

figure % Figure 3.6
tile=tiledlayout(3,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([3 1])
hold on
plot(I*100,mibl*100,'Color','k','LineStyle','-','LineWidth',0.5)
plot(I*100,mi*100,'Color','k','LineStyle','--','LineWidth',0.5)
%plot(RR,mov,'Color','k','LineStyle','--','LineWidth',0.5) 
box on
ytickformat('percentage');
xtickformat('percentage');
ylabel('Inflation rate ({\pi})','FontSize',20)
%title('Panel 1: Expected inflation in the estimated stationary distribution')
legend('Baseline','\sigma = 0.9','FontSize',15)
legend boxoff
xlim([0 7]);
xlabel('Interest rate ({\it{i}})','FontSize',20)
hold off

figure % Figure 3.6
tile=tiledlayout(3,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([3 1])
hold on
plot(I*100,cost,'Color','k','LineStyle','-','LineWidth',0.5)
box on
xtickformat('percentage');
ylabel('Expected Cost','FontSize',20)
%title('Panel 1: Expected cost in the estimated stationary distribution')
%legend boxoff
xlim([0 7]);
xlabel('Interest rate ({\it{i}})','FontSize',20)
hold off


%output for figure 3.6 from baseline
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\movbl mov -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\mvbl mv -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\pmvbl pmv -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\pRR pRR -ASCII -DOUBLE;

%output for figure 3.8 input to ASBconV2Fig38.m
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Rcon103 RR -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Gbar103 mgb -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\mabl103 ma -ASCII -DOUBLE;

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\TGB Tc -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\TBB Tu -ASCII -DOUBLE;

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\EDGB3 elsc -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\EDBB3 elsu -ASCII -DOUBLE;




