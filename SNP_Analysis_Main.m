%snp analysis!
%% outputting an approrpiate snp file from https://genome.ucsc.edu/cgi-bin/hgTables this will be manually fed into INDD
%sigh, what a pain, we have to jump in and out of IGB browser and INDD
outputRead=true;
if outputRead
rawReads=readtable('/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/RSBN1'); %this snp list comes from IGB
writetable(rawReads, '/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/RSBN1.xls')
end

%% Initializing reference variables
snpReference=readtable('/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/SNP_referenceSheet.xlsx'); %it's irritating to have to design the code this way, but INDD is just too fussy this sheet contains info like diagnosis date and gene status
geneName='FAF1';
%FAF1, ACTN2 , DNM1P46, KMT2C, KDM5B
%sposed to be negative "controls" RSBN1 ,**MIA, TLX1, Linc00624, COG5, 
[snpBasePath,snpOutputPath,geneSite, peakSite]=snpLookupTable(geneName);
snpTable=readtable(snpOutputPath); %this should come from the output .raw file

%/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/DNM1P46FullOutput.txt



idRef= snpReference.INDDID ; ageAtOnsetHold= snpReference.AgeatOnset; dx=snpReference.ClinicalPhenotype1; mutSum= snpReference.Mutation_Summary;
SNPstudy=snpReference.GWAStudy; snpOutputDate=snpReference.PANDoRA_Date;

% 1 for control 2 for ALS 3 for FTD 4 for ALS-FTD 
dxCodeHold=nan(1,length(dx));
for dxCount=1:length(dx)
    if contains(dx{dxCount}, 'ALS-FTD')
        dxCodeHold(dxCount)=4;
    elseif contains(dx{dxCount}, 'FTD')
        dxCodeHold(dxCount)=3;
    elseif contains(dx{dxCount}, 'ALS')
        dxCodeHold(dxCount)=2;
    elseif contains(dx{dxCount}, 'Normal')
        dxCodeHold(dxCount)=1;
    end
end

studyCodeHold=false(1,length(dx));

% for dd=1:length(snpOutputDate)
%     if ~isempty(snpOutputDate{dd})
%         % str2num(snpOutputDate{dd}(end-3:end))
% studyCodeHold(dd)= str2num(snpOutputDate{dd}(end-3:end)) <= 2015;
%     end
% end
% for dd=1:length(SNPstudy)
%      if contains(SNPstudy{dd}, 'Van Deerlin FTLD-TDP GWAS')
%          studyCodeHold(dd)=true;
%      % elseif contains(SNPstudy{dd}, 'ADCC GWAS Phase 9')
%      %    studyCodeHold(dd)=true;
%       % elseif contains(SNPstudy{dd}, 'Singleton')
%       %    studyCodeHold(dd)=true;
%      % elseif contains(SNPstudy{dd}, 'CNDR GSA 2019')
%      %           studyCodeHold(dd)=true;
%     % elseif contains(SNPstudy{dd},'NACC ALS FTD GDA 2022')
%     %                    studyCodeHold(dd)=true;
%   % elseif contains(SNPstudy{dd},'ADCC GWAS Phase 8')
%   %                      studyCodeHold(dd)=true;
%  % elseif contains(SNPstudy{dd},'ALS GSA 2019')
%  %          studyCodeHold(dd)=true;
%  % elseif contains(SNPstudy{dd},'AD-SYN GWAS')
%  %          studyCodeHold(dd)=true;
%  % elseif contains(SNPstudy{dd},'%ALS-GWAS')
%  %          studyCodeHold(dd)=true;
% % elseif isempty(SNPstudy{dd})
% %         studyCodeHold(dd)=true;
%     end
% end


geneCodeHold=false(1,length(mutSum)); %for the time being we'll only analyze ftd patients with prg and c9 and exclude ALS patients with SOD1 mutations
geneCodeHold(dxCodeHold==2 | dxCodeHold==4 |dxCodeHold==1)=true; %all ALS and control patients

for geneCount=1:length(geneCodeHold)
    if contains(mutSum{geneCount}, 'C9orf72') || contains(mutSum{geneCount}, 'GRN') %add in proganulin and c9 FTD cases
        geneCodeHold(geneCount)=true;
    elseif contains(mutSum{geneCount}, 'SOD1') || contains(mutSum{geneCount}, 'FUS') %remove SOD1/FUS ALS
                geneCodeHold(geneCount)=false; 
    end
end

%removing non relevant ids in the evaluated output
delIdx=false(1, height(snpTable));
% triple checking to make sure no duplicates



IDAtPlay=snpTable.INDDID;
for dd=1:length(IDAtPlay)
    if sum(idRef== IDAtPlay(dd) ) ~=1
        delIdx(dd)=true;
    end
end


snpTable(delIdx,:)=[];
IDAtPlay=snpTable.INDDID;

delIdx=false(1, height(snpTable));
% triple checking to make sure no duplicates

for dd=1:length(IDAtPlay)
if sum(IDAtPlay==IDAtPlay(dd))>1
    delIdx((IDAtPlay==IDAtPlay(dd)))=true;
    delIdx(find(IDAtPlay==IDAtPlay(dd),1))=false;
end
end
snpTable(delIdx,:)=[];
IDAtPlay=snpTable.INDDID;


%sigh ok last thing, reorder all the previously produced indexes obvi you
%can see I wrote this code out of order because I didn't anticipate how
%frustrating INDD would be to use
dxCode=nan(1,length(dx));
geneCode=false(1,length(mutSum)); %for the time being we'll only analyze ftd patients with prg and c9 and exclude ALS patients with SOD1 mutations
ageAtOnset=nan(1,length(idRef));
mutsumTrue=cell(size(mutSum));
studyCodeTrue=false(1,length(SNPstudy));
for dd=1:length(idRef)
newSpot=find(IDAtPlay==idRef(dd));
geneCode(newSpot)=geneCodeHold(dd);
dxCode(newSpot)=dxCodeHold(dd);
ageAtOnset(newSpot)=ageAtOnsetHold(dd);
mutSumTrue{newSpot}=mutSum{dd};
studyCodeTrue(newSpot)=studyCodeHold(dd);
end





%% analysis/plotting

snpDataALS=[]; snpDataFTD=[]; snpDataTot=[];

allCol=snpTable.Properties.VariableNames; %this will help us iterate through all snps

for varCount=1:length(allCol)
    if contains(allCol{varCount},'rs' )%if this is a SNP column
        snp2analyze=table2array(snpTable(:,varCount)); %output counts
        if iscell(snp2analyze(1))
            continue
        end
    else
        continue
    end

% if sum(isnan(snp2analyze(dxCode==1)))>=5 || sum(isnan(snp2analyze(dxCode>1)))>=50
%     continue
% end

% if sum(isnan(snp2analyze)) > (.1*length(snp2analyze))
%     continue
% end

%inputs to snp stat runner are gene code, the path patients, the control
%patients, and the snp at play
% FTD
[pvals,propDiff,~] = snpStatRunner(geneCode,dxCode==3 & ~studyCodeTrue,dxCode==1 & ~studyCodeTrue,snp2analyze,allCol{varCount});
snpDataFTD=[snpDataFTD; {allCol{varCount}}, {pvals}, {propDiff}];

% ALS
[pvals,propDiff,~] = snpStatRunner(geneCode,dxCode==2 & ~studyCodeTrue,dxCode==1 & ~studyCodeTrue,snp2analyze,allCol{varCount});
snpDataALS=[snpDataALS; {allCol{varCount}}, {pvals}, {propDiff}];

%all viable patients
[pvals,propDiff,~] = snpStatRunner(geneCode,dxCode>1 & ~studyCodeTrue,dxCode==1 & ~studyCodeTrue,snp2analyze,allCol{varCount});
snpDataTot=[snpDataTot; {allCol{varCount}}, {pvals}, {propDiff}];
end  





%% plotting

snpPlotter(snpDataALS,'ALS',geneName, snpBasePath, geneSite, peakSite)
snpPlotter(snpDataFTD,'FTD',geneName, snpBasePath, geneSite, peakSite)
snpPlotter(snpDataTot,'ALS and FTD',geneName, snpBasePath, geneSite, peakSite)
%%  manhattan Plot






%% play around space
% 
% 
% 
% 
% 
% snpOfInterest=table2array(snpTable(:,28))';
% 
% % figure
% % scatter(snpOfInterest(dxCode==2 &geneCode), ageAtOnset(dxCode==2 &geneCode))
% 
% 
% mutSumTrue(snpOfInterest==0 &dxCode==2 & geneCode)
% 
% nanmean(ageAtOnset(dxCode==2 & geneCode & snpOfInterest==1))
% 
% nanmean(ageAtOnset(dxCode==2 & geneCode & snpOfInterest==0))
% 
% 
% 
% containsMut=false(1,length(mutSum));
% containsC9=false(1,length(mutSum));
% for tt=1:length(mutsumTrue)
% containsMut(tt)= ~isempty(mutSumTrue{tt});
%     if contains(mutSumTrue{tt}, 'C9orf72')
%     containsC9(tt)=true;
%     end
% end
% 
% 
% nanmean(ageAtOnset(dxCode==2 & geneCode & snpOfInterest==1 &containsC9 ))
% nanmean(ageAtOnset(dxCode==2 & geneCode & snpOfInterest==1 &containsMut))
% 
% 
% 
% mutSumTrue(dxCode==2 & geneCode & snpOfInterest==0)
% 
% nanmean(ageAtOnset(dxCode==2 & geneCode & snpOfInterest==0))
