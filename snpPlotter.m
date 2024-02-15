function [] = snpPlotter(snpData,dxName,geneName, snpBasePath, geneSite, peakSite)


heteroCarryP=[];
homoCarryP=[];
allMinorAlleleP=[];
heteroCarryProp=[];
homoCarryProp=[];
allMinorAlleleProp=[];


% taking data back out of original format
for tt=1:size(snpData,1)
heteroCarryP=[heteroCarryP,snpData{tt,2}(1)];
homoCarryP=[homoCarryP,snpData{tt,2}(2)];
allMinorAlleleP=[allMinorAlleleP,snpData{tt,2}(3)];
heteroCarryProp=[heteroCarryProp,snpData{tt,3}(1)];
homoCarryProp=[homoCarryProp,snpData{tt,3}(2)];
allMinorAlleleProp=[allMinorAlleleProp,snpData{tt,3}(3)];
end


figure
subplot(2,2,1)
scatter(homoCarryProp,-log10(homoCarryP), 'b')
hold on
scatter(homoCarryProp(logical(fdr_bh(homoCarryP)) ) ,-log10(homoCarryP(logical(fdr_bh(homoCarryP)) )  ), 'm')
legend({'q>.05' ,'q<.05'})
title('Homozygous Minor Allele Carrier Ratio')
ylabel('-log10(p)')
xlabel(['Percent Difference in Carrier Frequency ', dxName, ' - Control' ])
subplot(2,2,2)
scatter(allMinorAlleleProp,-log10(allMinorAlleleP),'b')
hold on
try
scatter(allMinorAlleleProp(logical(fdr_bh(allMinorAlleleP)) ) ,-log10(allMinorAlleleP(logical(fdr_bh(allMinorAlleleP)) )  ), 'r')
legend({'q>.05' ,'q<.05'})
title('Minor Allele Carrier Ratio')
sgtitle(['Difference in ', geneName, ' Carrier Frequency between ', dxName, ' and Control pts'])
catch
    b=1;
end




% manhattan Plot

snpIDX=readtable(snpBasePath);
snpIDXName=snpIDX.name;
snpStart= geneSite(1);
snpEnd=geneSite(2);
snpLoc=nan(1,height(snpData));

for snpAtPlay=1:height(snpData)
nameAtPlay= snpData{snpAtPlay,1}; nameAtPlay=nameAtPlay(1: (find(nameAtPlay=='_',1)-1)) ; %we truncate off the lil annotation
    for snpIDXr=1:length(snpIDXName)
        if contains(snpIDXName{snpIDXr}, nameAtPlay)
snpLoc(snpAtPlay)= floor(nanmean( [ table2array(snpIDX(snpIDXr,2)), table2array(snpIDX(snpIDXr, 3))   ] ));
break 
        end
    end
end

subplot(2,2,3)
scatter((snpLoc- (snpStart)),-log10(homoCarryP),'k')
hold on
scatter(( snpLoc(logical(fdr_bh(homoCarryP))  ) - (snpStart) ) ,-log10(homoCarryP(logical(fdr_bh(homoCarryP)) )  ), 'm')
vline(0, 'k--');
vline(snpEnd-snpStart, 'k--');
    if ~isempty(peakSite)
    vline(peakSite-snpStart, 'g--')
    end



subplot(2,2,4)

scatter((snpLoc- (snpStart)),-log10(allMinorAlleleP),'b')
hold on
scatter(( snpLoc(logical(fdr_bh(allMinorAlleleP))  ) - (snpStart) ) ,-log10(allMinorAlleleP(logical(fdr_bh(allMinorAlleleP)) )  ), 'r')
vline(0, 'k--');
vline(snpEnd-snpStart, 'k--');
    if ~isempty(peakSite)
    vline(peakSite-snpStart, 'g--')
    end


