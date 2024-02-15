function [snpPath,snpOutputPath,geneSite, peakSite] = snpLookupTable(geneName)
%Simple lookup table to make things easy

if contains(geneName,'FAF1')
    snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/FAF1_initial_input';
    snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/FAF1_FullOutput.txt';
    geneSite= [50410865,50986429];
    peakSite=[50511711,50512211];
elseif contains(geneName,'RSBN1')
snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/RSBN1';
   snpOutputPath= '/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/RSBN1OutputFull.txt';
    geneSite=[113759229, 113815008];
    peakSite=[];
elseif contains(geneName, 'ACTN2')
    snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/actn2_full.xls';
    snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/actn2FullOutput.xlsx';
    geneSite= [236682591, 236768538];
    peakSite=[236714481	236714981];
elseif contains(geneName, 'KDM5B')
    snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/KDM5B_raw_input';
    snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/KDM5BFullOutput.txt';
    geneSite= [202720970, 202813684];
    peakSite=[202744340	202744840];
elseif contains(geneName, 'MIA')
    snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/MIA_raw_input.txt';
    snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/MIA_Output_Full.txt';
    geneSite= [222615600, 222670503];
    peakSite=[];
elseif contains(geneName,'TLX1')
    snpPath= '/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/TLX1_raw';
    snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/TLX1FullOutput.txt';
    geneSite=[101131299,101137789];
    peakSite=[];
elseif contains(geneName, 'COG5')
snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/COG5';
snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/COG5FullOutput.txt';
    geneSite=[107201743,107564514];
    peakSite=[];
elseif contains(geneName, 'KMT2C')
snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/KMT2C_initial_input';
snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/KMT2CFullOutput.txt';
    geneSite=[152119870,152451057];
    peakSite=[152261214,	152261714];
elseif contains(geneName, 'Linc00624')
snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/KMT2C_initial_input';
snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/KMT2CFullOutput.txt';
    geneSite=[147382191,147517875];
    peakSite=[];
elseif contains(geneName, 'DNM1P46')
    snpPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Inputs/DNM1P46_initial_input';
    snpOutputPath='/Users/pennerc/Documents/Data/Chen_Plotkin_Rotation_2024/SNP_analysis/INDD_Outputs/DNM1P46FullOutput.txt';
    geneSite=[99790155,99806927];
    peakSite=[99799701,	99800201];


end 