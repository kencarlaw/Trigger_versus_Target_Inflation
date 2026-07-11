clear
load mibl26.mat;
NN=29;
ic=0.0025;
lam=40;
lam2=1;

I=zeros(NN,1);cost=zeros(NN,1);

for i=1:NN
    I(i)=(i-1)*ic;
    cost(i)=lam2*I(i)+lam*mibl(i);%+lam2*(1+I(i))^1.5;
end

figure % Figure 3.6
tile=tiledlayout(3,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([3 1])
hold on
plot(I*100,cost,'Color','k','LineStyle','-','LineWidth',0.5)
box on
xtickformat('percentage');
ylabel('Cost of monetary polcy')
%title('Panel 1: Expected cost in the estimated stationary distribution')
%legend boxoff
xlim([0 7]);
%ylim([0.05 0.24])
xlabel('Interest rate ({\it{i}})')
hold off