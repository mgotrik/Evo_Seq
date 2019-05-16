%%Step 1: Import files and identify sequence/library info:

% Choose FASTQ Directory
%FastQDirectory = 'D:\Sequencing Data\BFDEFGH Reads\'; %Where FastQ files are
FastQDirectory = 'C:\Users\Michael\OneDrive\Work Documents\9. Data\Sequencing_Temp\';


%FastQFiles = dir(strcat(FastQDirectory,'*','.fastq')); %Names of the FastQ Files
%NumFiles=length(FastQFiles); %Delete?

FilterDir = [FastQDirectory, 'FilteredFastQ\']; %Where the files are saved

BlockSize = 30000; %Number of Sequences to analyze at a time
minqual = 20; %min quality to qualify as 'low quality'
maxpct = 5; %max percent of low quality bases

if 7~=exist(FilterDir,'dir') %
    mkdir(FilterDir); %Makes the FilteredFastQ folder if it doesnt exist
end
%%
[PDS] = ES_ImportFastQInfo(FastQDirectory);

%%
%Parent Pools is the [PoolIndex,IndexOfParent] for each pool that was
%sequenced
ParentPools = [ [1,0]  %File index of the current pool, and the parent index
                [2,1]   % Parent = 0 means that it is the parent pool
                [3,0]
                [4,0]
                [5,4]
                [6,5]
                [7,6]
                [8,6]
                [9,7]];
%LibIDs is the Library ID (from Libraries.xls) for each [n,0] index in the
%ParentPools
LibIDs = {[2],[4:6],[1,3]};

[PDS] = ES_getPDSLibInfo(PDS,ParentPools,LibIDs);

%%
%{
NumPools = size(PDS,2); %Number of pools that are being analyzed
%Insert Info on library
Pool1Files = 4:9;
Pool2Files = 1:2;
Pool3Files = 3;
%
for i = 1:NumPools
    if ismember(i,Pool1Files) %BFH Pools
        Libraries = {'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGACAGCAGTTCTAATACGACTCACTATAGGGATTGTGAGTGTCCGGCGCAAGCTCTAAAGAGTAACGGCCGACTCNNNNNNNNNNNNGAGATTTGGGCATTTGCAGTTCCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG', %BF-H-bc
            'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTTTTATGCTATAATTATTTCATGTAGTAAGGAGGTTGTATGGAAGACGTTCCTGGATCCGGGATTGTGAGTGTCCGGNNNNNNNNNNAACGTAACGGCCGACNNNNNNNNNNTTGGGCATTTGCAGTTCCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG'}; %BF-AR
    elseif ismember(i,Pool2Files) %BFD Pools
        Libraries = {upper('AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGACAGCAGTTCTAATACGACTCACTATAGGGATTGTGAGTGTCCGGnnnnnnctcacgccgcgaagaacgaagaacaagattgagtctgatggaaggacgaaACATGAGGATCACCCATGTttctgtctggacagactcaatcttgttcttcgggagcggcgtgagnnnnnnTTGGGCATTTGCAGTTCCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG')}; %BF-D-bc
    elseif ismember(i,Pool3Files) %BFEFG Pools
        Libraries = {'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGACAGCAGTTCTAATACGACTCACTATAGGGATTGTGAGTGTCCGGGCGATATGGATAGAAGGACGGGTCCCATATTGCGAGAAACGTGCGGAGTACATCGCGAAACTGGATCACAGGCGAGAAACAGACCGTGGCGGTAGTTNNNNNNNNNNAACTACCTTGAACGGTCTGGAGCCTGTGTGTTAGGAAACTAACTTTGTCCAGGAGATGGAAACGCACGGAGCAATATGGTTGAGTAGAGTGTGAGCTATCCTAAGTCGCTTGGGCATTTGCAGTTCCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG', %BF-E-bc
            'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGACAGCAGTTCTAATACGACTCACTATAGGGATTGTGAGTGTCCGGGCGATATGGATAGAAGGACGGGTCCCATATTGCGAGAAACGTGTGCAAGGCTCCCTGTGTTACAATGACGCTTCCACGTCGCGAAACTGNNNNNNNNNNNCAGGAGACGTGGAAGTAAGACAGTAACTTTGTAGCGGAAACGCTATCAGGGTAGCCACCACACGGAGCAATATGGTTGAGTAGAGTGTGAGCTATCCTAAGTCGCTTGGGCATTTGCAGTTCCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG', %BF-F-bc
            'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGACAGCAGTTCTAATACGACTCACTATAGGGATTGTGAGTGTCCGGGCGATATGGATAGAAGGACGGGTCCCATATTGCGAGAAACGTTGGCAGCCTACGAAGAACTGATCAACGGAGAAGTAGACTGCTACNNNNNNNNNNNGTAGAACTACAGTCTACTTCGAGAAACGTTGTGAGCGGAAACGCTCTTTGTCAGGGAGTAGGCAGCCAACGGAGCAATATGGTTGAGTAGAGTGTGAGCTATCCTAAGTCGCTTGGGCATTTGCAGTTCCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG' %BF-G-bc
            };
    end
    ForwardPrimer = 'CGACAGCAGTTCTAATACGACTCACTATAG'; %the library FP
    ReversePrimer = 'TGGAACTGCAAATGCCCAA'; %the library RP, not the RPcomp
    SequenceFP = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'; %Sequencing FP (Read 1 TrueSeq)
    SequenceRP = 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'; %Actual Sequencing RP (not ReverseComplement, Read 2 Trueseq)
    %Gets the library that is sequenced

    
    for j = 1:length(Libraries)
        TempLib = Libraries{j}; %Looks at library #i
        TempLib = strsplit(TempLib,SequenceFP);
        TempLib = strsplit(TempLib{2},seqrcomplement(SequenceRP));
        SeqLibrary = TempLib{1}; %Finds region that is sequenced
        
        SeqLength = PDS(i).Read1Len + PDS(i).Read2Len; %Gets total sequencing length
        OverlapLen = SeqLength - length(SeqLibrary); %Difference between length of library and how much was sequenced
        PDS(i).OverlapLen(j) = OverlapLen;
        
        R1MaxLen = min(PDS(i).Read1Len,length(SeqLibrary));
        PDS(i).SeqLibrary{j} = SeqLibrary;
        Library = SeqLibrary(strfind(SeqLibrary,ForwardPrimer):strfind(SeqLibrary,seqrcomplement(ReversePrimer))+length(ReversePrimer)-1);
        PDS(i).Library{j} = Library; %Finds region of interest (i.e. the library primers)
        PDS(i).LibraryLen(j) = length(Library);
        PDS(i).Seq_Read1{j} = SeqLibrary(1:R1MaxLen);
        OverlapStartR1 = max(1,PDS(i).Read1Len-OverlapLen+1); %Returns where in library the Read 2 reaches
        PDS(i).Seq_Overlap{j} = SeqLibrary(OverlapStartR1:R1MaxLen);
        PDS(i).Seq_Read2{j} = SeqLibrary(OverlapStartR1:end);
        PDS(i).Seq_Read2rc{j} = seqrcomplement(SeqLibrary(OverlapStartR1:end));        
    end
    
    %PoolDataStruct(i).Library = Libraries{1}
    PDS(i).BlockSize = BlockSize;
    PDS(i).FP = ForwardPrimer;
    PDS(i).RP = ReversePrimer;
    PDS(i).RPc = seqrcomplement(ReversePrimer);
    PDS(i).FPLen = length(ForwardPrimer);
    PDS(i).RPLen = length(ReversePrimer);
    RPc = seqrcomplement(ReversePrimer);
    PDS(i).maxFPScore = swalign(ForwardPrimer,ForwardPrimer,'GapOpen',8,'ExtendGap',10,'Alphabet','NT');
    PDS(i).maxRPcScore = swalign(RPc,RPc,'GapOpen',8,'ExtendGap',10,'Alphabet','NT');
    PDS(i).MutRP = 0;
    PDS(i).MutFP = 0;
    PDS(i).NoFP = 0;
    PDS(i).NoRP = 0;
    PDS(i).CorrBases = 0;
    PDS(i).AlignFail = 0;

    %Identifies whether to consider Read 1 or Read 2 based on presence of
    %FP and/or RP.
    for j = 1:length(Libraries)
        if isempty(strfind(PDS(i).Seq_Read1{j},ForwardPrimer))
            ReadFlag = 2; %Consider only Read 2
        elseif isempty(strfind(PDS(i).Seq_Read2{j},RPc))
            ReadFlag = 1; %Consider only Read 1
        elseif ~isempty(strfind(PDS(i).Seq_Read1{j},ForwardPrimer)) && ~isempty(strfind(PDS(i).Seq_Read2{j},RPc))
            ReadFlag = 3; %Library is covered in both Read1 and Read2
        else
            error('Logic failed, oh no!')
        end
        PDS(i).ReadFlag(j) = ReadFlag;
    end
end
%}

%% Just initiates parallel pool
parfor k=1
    disp('Parallel pool online');
end
%%
PDS = ES_GetAnalysisInfo(PDS,BlockSize,minqual,maxpct,FilterDir);
%{
if ~exist("PoolDataStruct(end).FilePathFilter") || 2~=exist(PDS(end).FilePathFilter{1, 1},'file') %Checks to see if the filtered fastq files already exist
%% Filters low quality reads and saves updated fastq files to FilterDir
    PDS(1).FilePathFilter = ''
    PDS(1).FailedQualityCheck = 0;
    PDS(1).FilteredEntries=0;
    PDS(1).NoSeqBlocks = 0; %NoSeqBlocks = floor(SeqNoStart(i)/SeqBlockSize);
    PDS(1).LastBlockSize = 0;
    tic
    parfor i=1:2
        [PDS(i).ADS.FilePathFilter,PDS(i).FilteredEntries,PDS(i).FailedQualityCheck] = seqfilter([{char(PDS(i).FilePath(1))},{char(PDS(i).FilePath(2))}],'PairedFiles',true,...
                                 'Method','MaxPercentLowQualityBases',...
                                 'Threshold',[minqual maxpct],...
                                 'OutputDir', FilterDir,...
                                 'UseParallel', true,...
                                 'WriteSingleton',false); %Does not write ones that have one read fail but one good
        PDS(i).FilteredEntries = PDS(i).FilteredEntries(1);
        PDS(i).FailedQualityCheck=PDS(i).FailedQualityCheck(1);
        PDS(i).NoSeqBlocks = ceil(PDS(i).FilteredEntries(1)/PDS(i).BlockSize(1)); %NoSeqBlocks = floor(SeqNoStart(i)/SeqBlockSize);
        PDS(i).LastBlockSize = mod(PDS(i).FilteredEntries(1),PDS(i).BlockSize(1));
        NumFailed = PDS(i).FailedQualityCheck(1)/PDS(i).FilteredEntries(1);
        disp(['Number of low quality reads removed (%): ',num2str(PDS(i).FailedQualityCheck(1)),' (',num2str(NumFailed*100),'). Pool:',PDS(i).PoolName{1}]);
    end
    toc %270 seconds for 6 pools, 45s for 1 more; 13s for 7 pools (PC) %24 seconds (PC) for BF-AR
    %9 pools read length 305 & 9 pools - 106.5s
    %save([FilterDir,'PoolDataStruct.mat'],'PoolDataStruct');
end
%%
%% Filters low quality reads and saves updated fastq files to FilterDir
    
    tic
    parfor i=1:10
        PDS(i).ADS.FilePathFilter = '';
        [PDS(i).ADS.FilePathFilter,PDS(i).ADS.FilteredEntries,PDS(i).ADS.FailedQualityCheck] = seqfilter([{char(PDS(i).ADS.FilePaths(1,1))},{char(PDS(i).ADS.FilePaths(1,2))}],'PairedFiles',true,...
                                 'Method','MaxPercentLowQualityBases',...
                                 'Threshold',[minqual maxpct],...
                                 'OutputDir', FilterDir,...
                                 'UseParallel', true,...
                                 'WriteSingleton',false); %Does not write ones that have one read fail but one good
        PDS(i).ADS.FilteredEntries = PDS(i).ADS.FilteredEntries(1);
        PDS(i).ADS.FailedQualityCheck=PDS(i).ADS.FailedQualityCheck(1);
        PDS(i).ADS.NoSeqBlocks = ceil(PDS(i).ADS.FilteredEntries(1)/PDS(i).ADS.BlockSize(1)); %NoSeqBlocks = floor(SeqNoStart(i)/SeqBlockSize);
        PDS(i).ADS.LastBlockSize = mod(PDS(i).ADS.FilteredEntries(1),PDS(i).ADS.BlockSize(1));
        NumFailed = PDS(i).ADS.FailedQualityCheck(1)/PDS(i).ADS.FilteredEntries(1);
        disp(['Number of low quality reads removed (%): ',num2str(PDS(i).ADS.FailedQualityCheck(1)),' (',num2str(NumFailed*100),'). Pool:',PDS(i).PoolName{1}]);
    end
    toc %270 seconds for 6 pools, 45s for 1 more; 13s for 7 pools (PC) %24 seconds (PC) for BF-AR
    %9 pools read length 305 & 9 pools - 106.5s
    %save([FilterDir,'PoolDataStruct.mat'],'PoolDataStruct');
%}
%%
save(['D:\Sequencing Data\BFDEFGH Reads\FilteredFastQ\PDS0514.mat'],'PDS')
%%
FastQDirectory = 'D:\Sequencing Data\BFDEFGH Reads\'; %Where FastQ files are
FilterDir = [FastQDirectory, 'FilteredFastQ\']; %Where the files are saved

%{
BlockSize = 30000; %Number of Sequences to analyze at a time
minqual = 20; %min quality to qualify as 'low quality'
maxpct = 5; %max percent of low quality bases
%}
%%
load(['D:\Sequencing Data\BFDEFGH Reads\FilteredFastQ\PDS0514.mat'],'PDS')
%%
parfor i=1:10
    fclose('all')
end
%%
load([FilterDir,'PoolDataStruct.mat'],'PoolDataStruct')
%% BEGIN COMPARING ACTUAL SEQUENCES
for i=1:3
    %%
    clearvars BIStart BIEnd BlockSize GapFailed NumSeqs Temp MisAligned
    PDS(i).ADS.AlignFail = 0;Align
    PDS(i).ADS.GappedReads = 0;
    PDS(i).ADS.FilePaths{3,1} = strrep(PDS(i).ADS.FilePaths{2,1},'_filtered','_filteredmerged');
    PDS(i).ADS.FilePathDescription{3,1} = {'Combines reads & aligns to lib'};
    %{
    
    parfor i = 1:9
    end
    
    %}
    tic
    %Determines which Read to consider (1, 2, or both=3)
    ReadFlag = min([PDS(i).LibStruct.ReadFlag]);  %DETERMINES WHETHER TO SEARCH FOR FP, RP, OR BOTH

    %ReadFlag=3;
    %Determines the block size:
    SeqBlockSize = PDS(i).ADS.BlockSize; %NUMBER TO COMPARE
    %TempBlockSize = SeqBlockSize;
    LastBlockSize = PDS(i).ADS.LastBlockSize;
    NumBlocks = double(PDS(i).ADS.NoSeqBlocks); %NUMBER OF SEQUENCING BLOCKS
    
    %Identifies where to save the merged sequences
    FilePaths = PDS(i).ADS.FilePaths(1,:);
    %FilePaths = PDS(i).ADS.FilePaths(2)
    %FilePaths = PDS(i).FilePathFilter; %WHERE COMPARED SEQUENCES WILL BE SAVED
    NewFilePath = strrep(strrep(PDS(i).ADS.FilePaths(2),'_R1','_R3'),'_filtered','_filtered_merged');
    
    %Variables that count how many sequences were missing things;
    NoPrimers=0; %COUNTER FOR SEQUENCES WITHOUT PRIMERS
    Correction=0; %COUNTER FOR NUMBER OF CORRECTIONS MADE (EACH BASE)
    warningID = 0;
    warnflag=0;
    
    %Defines regions of the Library
    %OverlapSeq = PDS(i).LibStruct.ReadOL{1}; %REGION THAT IS OVERLAPPED BETWEEN R1 AND R2
    FP = PDS(i).LibStruct(1).FP;
    RP = PDS(i).LibStruct(1).RP;
    RPc = seqrcomplement(PDS(i).LibStruct(1).RP);
    %%
    for j=1:NumBlocks
        MisAligned(j) = 0;
        GapFailed(j) = 0;
        Temp(j) = struct('Seqs',struct(),'SeqF',struct(),'SeqR',struct());%struct('Header',[],'Sequence',[],'Quality',[]);
        NumSeqs(j) = 0;
        if j==NumBlocks %IF WE ENCOUNTER THE LAST BLOCK, REDEFINE HOW MANY WE ARE LOOKING AT
            BlockSize(j)=LastBlockSize;
        else
            BlockSize(j) = SeqBlockSize;
        end
        BIStart(j) = (j-1)*SeqBlockSize+1; %Indexing where the block of sequences starts
        BIEnd(j) = (j-1)*SeqBlockSize+BlockSize(j); %indexing where the block of sequences end
    end
    %%
    %TempSeqF = cell2struct(cell(NumBlocks,BlockSize))
    %TempSeqR = cell(NumBlocks,BlockSize)
    %TempSeqs = cell(NumBlocks,BlockSize)
    tic
    parfor j=1:NumBlocks %GOES ONE BLOCK AT A TIME
        %%
    
        %{
        if j==NumBlocks %IF WE ENCOUNTER THE LAST BLOCK, REDEFINE HOW MANY WE ARE LOOKING AT
            TempBlockSize=LastBlockSize;
        end
        BIStart = (j-1)*SeqBlockSize+1; %Indexing where the block of sequences starts
        BIEnd = (j-1)*SeqBlockSize+TempBlockSize; %indexing where the block of sequences end
        
        %}
        if ReadFlag==3 %IF WE WANT TO LOOK AT BOTH READS
            %%
            %{a
            
            disp(['entering into clean:', num2str(j)]);
            Temp(j).SeqF = fastqread(char(FilePaths(1)),'blockread', [BIStart(j) BIEnd(j)]); %Loading Forward strand
            Temp(j).SeqR = fastqread(char(FilePaths(2)),'blockread', [BIStart(j) BIEnd(j)]); %Loading Reverse Strand
            % Aligns two sequences to library/each other
            Temp(j).SeqR = ES_QS_SeqRCompBlock(Temp(j).SeqR);

            [Temp(j).Seqs,AlignFail,GappedReads] = ES_B_CleanAlignBlocks(Temp(j).SeqF,Temp(j).SeqR,BlockSize(j));
            %}
            %[test1,test2,test3] = ES_B_CleanAlignBlocks(fastqread(char(FilePaths(1)),'blockread', [BIStart(j) BIEnd(j)]),ES_QS_SeqRCompBlock(fastqread(char(FilePaths(2)),'blockread', [BIStart(j) BIEnd(j)])),BlockSize(j));
            %[Temp(j).Seqs,MisAligned(j),GapFailed(j)] = ES_B_CleanAlignBlocks(fastqread(char(FilePaths(1)),'blockread', [BIStart(j) BIEnd(j)]),ES_QS_SeqRCompBlock(fastqread(char(FilePaths(2)),'blockread', [BIStart(j) BIEnd(j)])),BlockSize(j));
            Temp(j).SeqF = ''
            Temp(j).SeqR = ''
            %clearvars test
            disp(['finished cleaning:', num2str(j)]);
            NumSeqs(j) = length(Temp(j).Seqs);
            %MisAligned(j) = AlignFail;
            %GapFailed(j) = GappedReads;
            %PDS(i).ADS.AlignFail = PDS(i).ADS.AlignFail
            %[TempSeq,CountFail,CountSeqCons] = PairReads05062019(TempSeqF,TempSeqR,PDS(i));
        %{    
        elseif ReadFlag==2 %IF WE WANT TO LOOK AT READ 2 ONLY
            %{
            TempSeqs = fastqread(char(FilePaths(2)),'blockread', [BIStart BIEnd]); %Loading Reverse Strand
            for k = 1:TempBlockSize %HAVE TO REVERSE THE SEQUENCES IN READ 2
                TempSeqs(k).Sequence = seqrcomplement(TempSeqs(k).Sequence);
                TempSeqs(k).Quality = fliplr(TempSeqs(k).Quality);
            end
            %}
            %%
            TempSeqs(j,:) = fastqread(char(FilePaths(2)),'blockread', [BIStart BIEnd]); %Loading Reverse Strand
            TempSeqs(j,:) = ES_QS_SeqRCompBlock(TempSeqs(j,:));
            %GlobalFPStart = strfind(PoolDataStruct(i).Seq_Read2,FP);
            %GlobalRPStart = strfind(seqrcomplement(PoolDataStruct(i).Seq_Read2),RP);   
        elseif ReadFlag==1 %IF WE WANT TO LOOK AT READ 1 ONLY
            TempSeqs(j,:) = fastqread(char(FilePaths(1)),'blockread', [BIStart BIEnd]); %Loading Forward strand
            %GlobalFPStart = strfind(PoolDataStruct(i).Seq_Read1,FP);
            %GlobalRPStart = strfind(seqrcomplement(PoolDataStruct(i).Seq_Read1),RP);
            %TempSeqR = fastqread(char(FilePaths(2)),'blockread', [BIStart BIEnd]); %Loading Reverse Strand
            %}
        end
        %% Will want to initialize this and have it look at all libraries
        %[TempSeqs,PDS(i)] = SeqScript_201905_CleanSequenceBlock(TempSeqs, PDS(i),TempBlockSize);
        %%
        %{
        if j==1
            Sequences=TempSeqs;
        else
            Sequences=horzcat(Sequences, TempSeqs);
        end
        %}

        warningID = warning('query','last');
        warning('off',warningID.identifier);
        
        %Temp(j).Seqs = [];
        %PoolDataStruc(i).FailedAlignment = AlignFail;
        disp(['pool' int2str(i) 'block' int2str(j)])
        
    end
    
    
    
    FinalStruct = Temp(1).Seqs;
    for j=2:NumBlocks
        FinalStruct = [FinalStruct; Temp(j).Seqs];
        Temp(j).Seqs = '';
    end
    fastqwrite(PDS(i).ADS.FilePaths{3,1},FinalStruct);
    
    clearvars FinalStruct
    PDS(i).ADS.AlignFailed = sum([MisAligned]);
    PDS(i).ADS.GapFailed = sum([GapFailed]);
    PDS(i).ADS.NumSeqsAfterAlign = sum([NumSeqs]);
    
    
    toc
    
    %%
    %{
    %PoolDataStruct(i).NoSeqAlignFail=AlignFail;
    %PoolDataStruct(i).NoSeqPrimersMissing=NoFP + NoRP;
    %PoolDataStruct(i).CorrBases=Correction;
    MergeFilename = strrep(PDS(i).FilePathFilter,'.fastq','_MergeAlign.fastq');
    PDS(i).MergeFilename=char(MergeFilename(1));
    fastqwrite(char(PDS(i).MergeFilename),Sequences);
    clearvars Sequences
    toc
    %}
end
save([FilterDir,'PDS_filtered_merged0514.mat'],'PDS');
%%
if exist('PoolDataStruct','var')==0 %Loads the PoolDataStruc if it is not loaded
    FilterDir = 'D:\Sequencing Data\BFAR1 Reads\FilteredFastQ\';
    load([FilterDir,'PoolDatastruc_filtered_merged.mat'])
end
%%
for i=fliplr(1:3)
    %Loads in Sequence file from FastQ. 
    %Removes fields except Sequences and adds 'marked' and 'count'
    Sequences = fastqread(PDS(i).ADS.FilePaths{3});
    Sequences = rmfield(Sequences,'Header');
    Sequences = rmfield(Sequences,'Quality');
    clearvars SeqMap KEYS COUNTS
    flag1=0;    %Set to 0 for the first run of a new pool
    StructSize = length(fieldnames(Sequences));
    Seqcell = reshape(struct2cell(Sequences),StructSize,[]);
    Seqcell = vertcat(Seqcell,num2cell(ones(1,length(Sequences))));
    Sequences = cell2struct(Seqcell,[fieldnames(Sequences);'marked'],1);
    
    StructSize = length(fieldnames(Sequences));
    Seqcell = reshape(struct2cell(Sequences),StructSize,[]);
    Seqcell = vertcat(Seqcell,num2cell(ones(1,length(Sequences))));
    Sequences = cell2struct(Seqcell,[fieldnames(Sequences);'count'],1);
    clearvars Seqcell;
    %Initializes the container size
    Div_1 = 3000;
    tic
    %Loops until all sequences have been counted
    countvar=0;
    while ge(length(Sequences),1)==1
        
        LenSeq = length(Sequences);
        countvar=countvar+1;
        PDS(i).LenSeq(countvar)=LenSeq;
        %Initializes the container with first unmarked sequence (all unmarked at this time)
        
        disp(['start ', num2str(length(Sequences))]);
        SeqMap = containers.Map(Sequences(1).Sequence,Sequences(1).count);
        Sequences(1).marked = str2num('0'); %marks first sequence
        j=2; %begins counting at i=2;
        
        % Creates a container of size Div_1 (e.g. 3000) sequences
        while le(length(SeqMap),Div_1) && le(j,LenSeq)
            %Checking i < # of sequences (e.g. we can still fill container)
            %{
            if gt(i,LenSeq) 
                break
            end
            %}
            %Checking if sequence is in container
            SeqTemp=Sequences(j).Sequence; 
            if isKey(SeqMap,SeqTemp)==0 %if it does not exist in the map
                TempCount = Sequences(j).count;
                SeqMap=[SeqMap; containers.Map(SeqTemp,TempCount)];
            else %if it does exist in the map
                TempCount = Sequences(j).count + SeqMap(SeqTemp); %counts occurences of sequence
                SeqMap=[SeqMap; containers.Map(SeqTemp,TempCount)]; %updates the occurence count in map
            end
            Sequences(j).marked = str2num('0');
            j=j+1;
            %{
            if isequal(i,LenSeq)
                break
            end
            %}
        end
        %j=i;
        % Begins comparing all other sequences to the first block of 3000
        while lt(j,length(Sequences))
            if mod(j,50000)==0
                disp([num2str(j),'/',num2str(length(Sequences)),' ',num2str(toc),'s']);
                tic
            end
            SeqTemp=Sequences(j).Sequence; %sequence thats being compared
            if isKey(SeqMap,SeqTemp)==0 %if it does not exist in the map
                j=j+1;
            else %if it does exist in the map, adjusts count in the map, marks sequence
                TempCount = Sequences(j).count + SeqMap(SeqTemp); %counts occurences of sequence
                SeqMap=[SeqMap; containers.Map(SeqTemp,TempCount)]; %updates the occurence count in map
                Sequences(j).marked = str2num('0');
                j=j+1;
            end
            %runs until all sequences in container have been mapped
        end
        if flag1==0 %if first run through, marks down the keys/values
            KEYS=SeqMap.keys;
            COUNTS = SeqMap.values;
            flag1=1;
            clearvars SeqMap
        else
            KEYS = horzcat(KEYS,SeqMap.keys); %cell containing all sequences made to length of block
            COUNTS = horzcat(COUNTS,SeqMap.values); %cell containing all copy numbers made to length of block
            clearvars SeqMap
        end
        %Sorts sequences based on 'marked' as above
        FieldNo = length(fieldnames(Sequences));
        Seqcell = reshape(struct2cell(Sequences),FieldNo,[]);
        Seqcell = Seqcell(1:FieldNo,:);
        CountFieldNames = fieldnames(Sequences);
        for k=1:length(CountFieldNames);
            if strcmp(CountFieldNames(k),'marked')==1
                CountField = k;
            end
        end
        Seqcell = flipud(sortrows(Seqcell',CountField))';
        Sequences = cell2struct(Seqcell,CountFieldNames,1);
        for k=1:length(Sequences)
            if Sequences(k).marked==0
                Sequences(k:end)='';
                break
            end
        end
        
        % check to make sure no sequences are being lost
        disp(['end ', num2str(length(Sequences) + sum(cell2mat(COUNTS)))]);
    end
%   Assigns Sequences to structure with new field 'count'.  Current length
%       of sequences is zero.
    for j=1:length(COUNTS) %reassigns sequences to a copy number
        Sequences(j).Sequence = char(KEYS(j)); %assigns actual sequence to map sequence
        Sequences(j).count = cell2mat(COUNTS(j)); %assigns actual copy number to map number
    end
%   Sorts sequences based on copy number in descending order
    Sequences = rmfield(Sequences,'marked');
    FieldNo = length(fieldnames(Sequences));
    Seqcell = reshape(struct2cell(Sequences),FieldNo,[]);
    Seqcell = Seqcell(1:FieldNo,:);
    CountFieldNames = fieldnames(Sequences);
    for j=1:length(CountFieldNames);
        if strcmp(CountFieldNames(j),'count')==1
            CountField = j;
        end
    end
    Seqcell = flipud(sortrows(Seqcell',CountField))';
    Sequences = cell2struct(Seqcell,CountFieldNames,1);
%   Saves Sequences variable to a .mat file
    disp(['Ending Sum is ',num2str(sum([Sequences(:).count]))]);
    PDS(i).ADS.FilePaths{4} = char(strcat(PDS(i).FileDir,PDS(i).PoolName,'_counted.mat'));
    PDS(i).VarName{1} = ['Seqs_',PDS(i).PoolName{1}];

    S = struct;
    S.(PDS(i).VarName{1}) = Sequences;
    save(PDS(i).ADS.FilePaths{4},'-struct', 'S');        %%ADJUST%%
    clearvars -except PDS FilterDir
    %clear all
end
%toc
save([FilterDir,'PDS_filtered_merged_counted.mat'],'PDS');
%% Calculates the fitness data.  Requires PoolDatastruct_filtered_merged_counted.mat. Assigns new .mat file to PoolDataStruct(i).SavedMatFileFitness
clearvars -except PDS
    CurrThreshold = 100;
    PrevThreshold = 10;
[SeqStruct,PDS] = ArrangePoolsForSequenceAnalysis(PDS,CurrThreshold,PrevThreshold);
save([FilterDir,'PDS_filtered_merged_counted_fitness.mat'],'PDS');
%% Makes a shorter fitness file for easier analysis
for i=4:9
    S = struct;
    CurrFieldName = PDS(i).PoolName{1}
    S.(PDS(i).PoolName{1}) = load(PDS(i).SavedMatFileFitness,CurrFieldName);
    S.(PDS(i).PoolName{1}) = S.(PDS(i).PoolName{1}).(PDS(i).PoolName{1});
    PDS(i).SavedMatFileFitnessShort=strrep(PDS(i).SavedMatFileFitness,'_fitness','_fitnessshort');
    
    for j =1:length(S.(PDS(i).PoolName{1}))
        if S.(PDS(i).PoolName{1})(j).count<10
            S.(PDS(i).PoolName{1})(j:end) = '';
            break
        end
    end
    save(PDS(i).SavedMatFileFitnessShort,'-struct', 'S');
    disp(num2str(i))
end
%%
    S= struct
for i=4:9
    CurrFieldName = PDS(i).PoolName{1}
    S.(PDS(i).PoolName{1}) = load(PDS(i).SavedMatFileFitnessShort,CurrFieldName);
    S.(PDS(i).PoolName{1}) = S.(PDS(i).PoolName{1}).(PDS(i).PoolName{1});
end
%%
count=0
sum(double([S.BFAR3ep(:).fitness]))
for i = 1:3341
    if isempty(S.BFAR3ep(i).fitness)
        count=count+1;
    end
end

%%

load('D:\Sequencing Data\BFAR1 Reads\FilteredFastQ\PoolDatastruct_filtered_merged_counted_fitness.mat')
%{
for i = 1:7
   
    FileName = char(strcat(PoolDataStruct(i).FileDir,PoolDataStruct(i).PoolName,'_counted.mat'))
    VarName = ['Seqs_',PoolDataStruct(i).PoolName{1}];
    PoolDataStruct(i).SavedMatFile = FileName;
    PoolDataStruct(i).SavedMatFileStructName = VarName;
end
%}
%% Groups sequences into Families based on similarity
%adds 'marked', 'Family' to structure field
SequenceTemp = Seqs_BFDR2;                                                    %%ADJUST%%
SeqSize = length(fieldnames(SequenceTemp));
Seqcell = reshape(struct2cell(SequenceTemp),SeqSize,[]);
Seqcell = vertcat(Seqcell,num2cell(ones(1,length(SequenceTemp))));
SequenceTemp = cell2struct(Seqcell,[fieldnames(SequenceTemp);'marked'],1);

SeqSize = length(fieldnames(SequenceTemp));
Seqcell = reshape(struct2cell(SequenceTemp),SeqSize,[]);
Seqcell = vertcat(Seqcell,num2cell(zeros(1,length(SequenceTemp))));
SequenceTemp = cell2struct(Seqcell,[fieldnames(SequenceTemp);'Family'],1);
clearvars Seqcell
l=1; %Sequence for comparison counter
i=1; %New sequence counter
k=1; %Family counter
%%
flag1=0;
while i<250%length(SequenceTemp)
    %   If sequence hasn't been seen yet, makes a new family out of it
    if SequenceTemp(i).marked ==1
        SequenceTemp(i).Family = k;
        FamilyTemp(k,l) = SequenceTemp(i);
        SequenceTemp(i).marked = str2num('0');
        %   Compares all other sequences to it, groups them into families
        %       if they compare
        for j=i+1:250%length(SequenceTemp)
            if SequenceTemp(j).marked ==1
                maxscore = nwalign(SequenceTemp(i),SequenceTemp(i));
                if gt(nwalign(SequenceTemp(i),SequenceTemp(j)),maxscore*0.9)==1
                    SequenceTemp(j).Family = k;
                    l=l+1;
                    FamilyTemp(k,l)=SequenceTemp(j);
                    SequenceTemp(j).marked = str2num('0');
                end
            end
        end
        %  Each cell of Family contains all members of the family, and all
        %  structural fitness variables
        NonEmpty = nnz(~cellfun('isempty',{FamilyTemp(k,:).Sequence}));
        Family(k)={FamilyTemp(k,1:NonEmpty)};
        k=k+1;
        l=1;
        disp(['Onto Family ',num2str(k)]);
    end
    i=i+1;
end
SequenceTemp = rmfield(SequenceTemp,'marked');
%   Creates a Fitness plot for the family to mark all single and double
%   point mutations vs. determined fitness

%%

%%
for i=2:size(Family,2) %steps through each member of family
    NonEmpty = nnz(~cellfun('isempty',{Family{1,i}(:).Sequence})); %# of sequences in family i
    if NonEmpty > 2 %if enough sequences to define a family
        ma = multialign({Family{1,i}(1:NonEmpty).Sequence},'terminalGapAdjust',true);
        FamilyConsensus{i}=cellstr(seqconsensus(ma,'gaps','all')); %sequence consensus of family
        ma = multialign({Family{1,i}(1:NonEmpty).Sequence,char(FamilyConsensus{i})},'terminalGapAdjust',true); %aligns all sequences with consensus
        FamConsTemp = ma(end,:); %
        FamConsTemp = strrep(FamConsTemp,'-','');
    %   Aligns ATCG up for each base, removes the current base from axis
        PlotMin = 'ATCG-'; %add/remove - to/to not include deletions
        PlotAxis = PlotMin; %gives the plot its axis
        AxisMin = length(PlotAxis);
        PlotAxis(1:AxisMin) = strrep(PlotAxis(1:AxisMin),FamConsTemp(1),'x');
        for j = 2:length(FamConsTemp)
            PlotAxis = [PlotAxis, PlotMin]; %concatenates ATCG for every base
            PlotAxis((j-1)*AxisMin+1:(j-1)*AxisMin+AxisMin) = strrep(PlotAxis((j-1)*AxisMin+1:(j-1)*AxisMin+AxisMin),FamConsTemp(j),'x');
        end
        PlotAxis = strrep(PlotAxis,'x','');
        AxisMin=AxisMin-1;
    %   Zeroes all parts of the matrix
        for k = 1:length(FamilyConsensus{1,i}{1})*AxisMin
            for l = 1:length(FamilyConsensus{1,i}{1})*AxisMin
                ScoringMatrix{i}(k,l) = 0;
            end
            ScoringMatrix{i}(k,k) = 0;
        end
    %   Assigns each single/double point mutation of matrix to a fitness
    %   level
        FamilyCNTot = sum([Family{1,i}(:).count]);
        for j=1:NonEmpty
            %Aligns each member of the family to consensus
            LocAlign = localalign(ma(j,:),ma(end,:),'GapOpen',10); %Aligns each to consensus
            MutLoc = regexp(LocAlign.Alignment{1}(2,:),'[^|]'); %Finds mutation locations
            MutLoc = setxor(regexp(LocAlign.Alignment{1}(3,:),'[-]'),MutLoc); %Ignores spaces
            Fitness = Family{1,i}(j).Fitness;
            FitnessWeighted = Fitness * (1 + Family{1,i}(j).count/FamilyCNTot)^2;
            %Counts if it is a single/double point mutation and assigns the
            %x/y value.  
            if length(MutLoc)==1
                xy = strfind(PlotAxis((MutLoc-1)*AxisMin+1:(MutLoc-1)*AxisMin+AxisMin),LocAlign.Alignment{1}(1,MutLoc));
                xy = xy + MutLoc*AxisMin;
                disp([num2str(j),' ',num2str(MutLoc),' ',num2str(xy),' ', LocAlign.Alignment{1}(1,MutLoc), ' ',LocAlign.Alignment{1}(3,MutLoc)]);
                ScoringMatrix{i}(xy,xy)=FitnessWeighted;
            elseif length(MutLoc)==2
                x =  strfind(PlotAxis((MutLoc(1)-1)*AxisMin+1:(MutLoc(1)-1)*AxisMin+AxisMin),LocAlign.Alignment{1}(1,MutLoc(1)));
                x = x + MutLoc(1)*AxisMin;
                y =  strfind(PlotAxis((MutLoc(2)-1)*AxisMin+1:(MutLoc(2)-1)*AxisMin+AxisMin),LocAlign.Alignment{1}(1,MutLoc(2)));
                y = y + MutLoc(2)*AxisMin;
                ScoringMatrix{i}(y,x)=FitnessWeighted;
            end
        end
        %
    
        
    end
end
%% Plots all family fitness plots
i=1
while i<=length(ScoringMatrix)
    if isempty(ScoringMatrix{i})==0
        figure
        imagesc(ScoringMatrix{i}(:,:,1))
        title(['Fitness Scores for Mutations of Family ',num2str(i)])
        xlabel('Sequence')
        ylabel('Sequence')
        colorbar
    end
    i=i+1
end
%% Populates Family variable with all properties
FamilyConsensus{2,i} = '';
for i = 1:length(Family)
    if length(Family{1,i})==length(FamilyConsensus)
        break
    end
    if i<=length(Family) && length(Family{1,i})>2
        Family{7,i} = FamilyConsensus{1,i};
        for j = 1:length(Family{1,i})
            if strcmp(Family{1,i}(j).Sequence, char(strrep(FamilyConsensus{1,i},'-','')))==1
                FamilyConsensus{2,i} = [Family{1,i}(j)];
                break
            end
        end
        if isempty(FamilyConsensus{1,i})==0 && isempty(FamilyConsensus{2,i})==0 
            Family{8,i} = FamilyConsensus{2,i};
            Family{9,i} = FamilyConsensus{2,i}.Fitness;
    end
    end
    Family{2,i} = sum([Family{1,i}(:).count]);
    Family{3,i} = sum([Family{1,i}(:).count])/sum([SequenceTemp.count]);
    Family{4,i} = max([Family{1,i}(:).Fitness]);
    Family{5,i} = median([Family{1,i}(:).Fitness]);
    Family{6,i} = min([Family{1,i}(:).Fitness]);
end
%%
MG3R4Fit = SequenceTemp;
MG3R4Family = Family;
%MG2R51ScoreMat = ScoringMatrix;
TopNo = 5;
FitnessSeqs = cell(TopNo,1);
CopyNoSeqs = cell(TopNo,1);
FamConsens = cell(TopNo,1);
for i = 1:TopNo
    FitnessSeqs(i) = {Family{1,i}(1).Sequence};
    CopyNoSeqs(i) = {SequencesMG3R4(i).Sequence};
    FamConsens(i)={strrep(Family{7,i},'-','')};
end

if ismember('OrderSeqs',who) == 0
    OrderSeqs = union(union(char(FamConsens{:}),CopyNoSeqs),FitnessSeqs)
else
    OrderSeqs = [OrderSeqs; union(union(char(FamConsens{:}),CopyNoSeqs),FitnessSeqs)]
end

