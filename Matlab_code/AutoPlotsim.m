clear

load sim26.mat;
load simsub1.mat;
load simsub2.mat;

figure %figure 3.1
tile=tiledlayout(3,2);
tile.Padding='none';
tile.TileSpacing='tight';
set(0,'DefaultAxesTitleFontWeight','normal');
nexttile([2 2])
autocorr(insim)
nexttile ([1 1])
autocorr(insub1)
nexttile([1 1])
autocorr(insub2)

[acfall]=autocorr(insim,NumLags=60);
[acf1]=autocorr(insub1,NumLags=60);
[acf2]=autocorr(insub2,NumLags=12);

dacfall=zeros(60,1);dacf1=zeros(60,1);
for i=1:60
    if i>1
        dacfall(i)=(acfall(i-1)-acfall(i))/acfall(i-1);
        dacf1(i)=(acf1(i-1)-acf1(i))/acf1(i-1);
    end
end
dacf2=zeros(12,1);
for i=1:12
    if i>1
        dacf2(i)=(acf2(i-1)-acf2(i))/acf2(i-1);
    end
end

macfall=mean(dacfall(1:40));
macf1=mean(dacf1(1:40));
macf2=mean(dacf2);
sacfall=std(dacfall(1:40));
sacf1=std(dacf1(1:40));
sacf2=std(dacf2);