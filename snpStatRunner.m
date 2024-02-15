function [pvals,propDiff,chiSquare] = snpStatRunner(geneCode,dxCode,contCode,snp2analyze,nameAtPlay)


  snp2analyze=snp2analyze(randperm(length(snp2analyze)));


%we count up heterozygotes homozygotes and all carriers of minor allele for
%path of interest and control
pathCounts=[sum(snp2analyze(dxCode & geneCode)==0),sum(snp2analyze(dxCode & geneCode)==1), sum(snp2analyze(dxCode & geneCode)==2) ];
contCounts= [sum(snp2analyze(contCode)==0),sum(snp2analyze(contCode)==1), sum(snp2analyze(contCode)==2) ];


%just outputs prop differences and runs a chi square test


propDiff(1)= ( sum(pathCounts(2)) /sum(pathCounts) )   - ( sum(contCounts(2)) /sum(contCounts) ) ;
propDiff(2)= ( sum(pathCounts(3)) /sum(pathCounts) )   - ( sum(contCounts(3)) /sum(contCounts) ) ;
propDiff(3)= ( sum( pathCounts(2:3)) /sum(pathCounts) )   - ( sum(contCounts(2:3)) /sum(contCounts) );

[~,pvals(1),chiSquare(1)]=  prop_test([pathCounts(2), contCounts(2)], [sum(pathCounts), sum(contCounts)],false);
[~,pvals(2),chiSquare(2)]=  prop_test([pathCounts(3), contCounts(3)], [sum(pathCounts), sum(contCounts)],false);
[~,pvals(3),chiSquare(3)]=prop_test([pathCounts(1), contCounts(1)], [sum(pathCounts), sum(contCounts)],false);


% sampleSize= sum(contCounts);
% chiHolderHomo=zeros(1,1000);
% ChiHolderHet=zeros(1,1000);
% ChiHolderTot=zeros(1,1000);
% trueVals=snp2analyze(dxCode & geneCode);
% for tt=1:200
% valPerm=trueVals(randperm(length(trueVals))); valPerm=valPerm(1:sampleSize);
% pathCounts=[sum(valPerm==0),sum(valPerm==1), sum(valPerm==2) ];
% [~,~,chiHolderHomo(tt)]=  prop_test([pathCounts(2), contCounts(2)], [sum(pathCounts), sum(contCounts)],false);
% [~,~,chiHolderHet(tt)]=  prop_test([pathCounts(3), contCounts(3)], [sum(pathCounts), sum(contCounts)],false);
% [~,~,chiHolderTot(tt)]=prop_test([(pathCounts(2)+pathCounts(3)), (contCounts(2)+contCounts(3))], [sum(pathCounts), sum(contCounts)],false);
% end
% 
% chiHolderHomo=nanmean(chiHolderHomo);
% chiHolderHet=nanmean(chiHolderHet);
% chiHolderTot=nanmean(chiHolderTot);
% 
% 
%         pvals(1) = 1 - chi2cdf(chiHolderHomo,1);
%         pvals(2) = 1 - chi2cdf(chiHolderHet,1);
%         pvals(3) = 1 - chi2cdf(chiHolderTot,1);



% 
%     pathCounts
%     contCounts
% % 
% if pvals(3)<.001
%     b=1;
%     pvals(2)
%     pathCounts
%     contCounts
%     nameAtPlay
% end


if isnan(pvals(2))
    pvals(2)=1; %this occurs when both have zero homozygotes
end

if isnan(pvals(3))
    pvals(3)=1; %this occurs when both have zero homozygotes or heterozygotes... should prob just remove these tbh... 
end


chiSquare=nan;
% [pvals,propDiff,chiSquare]=snpStatOutput(pathCounts,contCounts);






end