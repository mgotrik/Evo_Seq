function [PoolDataStruc] = ES_ImportFastQInfo(FastQDirectory)
%%
    %Edited 5/13/2019
    
    %%Identify FastQ files
    FastQFiles = dir(strcat(FastQDirectory,'*','.fastq')); %Names of the FastQ Files

    
    NumFiles = length(FastQFiles);
    for i=1:length(FastQFiles)
        SplitNames(i,:) = strsplit(FastQFiles(i).name,'_'); %Splits names into components
    end
    
    NameRoundIndex = 1;  %Index in fastqname that has Round/pool name
    NamePoolIDIndex = 2; %Index of S_ says which pool it is;
    NameReadIndex = 4; % index in fastqname that has the read
    
    %Identifies info from SplitNames
    PoolNames = unique(SplitNames(:,1)); %Identifies names of pools (e.g. BDAR0)
    NumPools = length(PoolNames); %Total number of pools
    ReadNames= unique(SplitNames(:,NameReadIndex));
    UniqueIDs = unique(SplitNames(:,2));
    for i=1:length(UniqueIDs)
        PoolID(i) = str2num(strrep(UniqueIDs{i},'S',''));
    end
    
    NumReads = length(ReadNames); %Number of Sequencing Reads (e.g. 1, 2)
    FilePaths = cell(NumPools,NumReads);

    %Indexing File paths into an array. Rows are Rounds, Columns are Reads
   
    %% Assembling information for each sequenced pool
    for j=1:NumPools %Counting along each round. e.g. Round1...Round2...
        i = PoolID(j);
        PoolDataStruc(i).ID = i;
        PoolDataStruc(i).PoolName = PoolNames(j);
        PoolDataStruc(i).FileDir=FastQDirectory;
        %Identify FileNames with this Pool Name
        Namefind = find(strcmp(SplitNames(:,NameRoundIndex),PoolNames(j)));
        
        %Identify the intersect of FileNames also with R1 and R2
        R1find = find(strcmp(SplitNames(:,NameReadIndex),'R1'));
        R2find = find(strcmp(SplitNames(:,NameReadIndex),'R2'));
        R1index = intersect(R1find,Namefind);
        R2index = intersect(R2find,Namefind);
        
        %{
        for j=1:NumReads %Counting along each Read. e.g. Read1...Read2...
            Readfind = find(strcmp(SplitNames(Namefind,4), ReadNames(j))); %index j for Read
            %FilePaths{i,j} = strcat(FastQDirectory,FastQFiles(Namefind(Readfind)).name);
            FastQInfoStruct(i,j) = fastqinfo(char(FilePaths(i,j))); %Imports info on FastQ Files
        end
        %}
        %Gets FileName for R1 and R2 for this pool
        FilePaths{i,1} = strcat(FastQDirectory,FastQFiles(R1index).name);
        FilePaths{i,2} = strcat(FastQDirectory,FastQFiles(R2index).name);
        %Gets FastQInfo 
        FastQInfoStruct(i,1) =  fastqinfo(char(FilePaths(i,1)));
        FastQInfoStruct(i,2) =  fastqinfo(char(FilePaths(i,2)));
        Read1 = fastqread(FilePaths{i,1},'blockread', [1000 1100]);
        Read2 = fastqread(FilePaths{i,2},'blockread', [1000 1100]);
        
        PoolDataStruc(i).LibStruct = struct;
        PoolDataStruc(i).FastqPath=FilePaths(i,:)'; %puts paths into a single structure
        PoolDataStruc(i).NumSeqs=FastQInfoStruct(i,:).NumberOfEntries; %Each one has same number of entries
        PoolDataStruc(i).ReadLen(1) = length(Read1(1).Sequence);
        PoolDataStruc(i).ReadLen(2) = length(Read2(1).Sequence);
        PoolDataStruc(i).FastQInfo = FastQInfoStruct(i,:);
        PoolDataStruc(i).ExampleSeqs(1:100,1) = {Read1(1:100).Sequence};
        PoolDataStruc(i).ExampleSeqs(1:100,2) = {Read2(1:100).Sequence};
        PoolDataStruc(i).ConsensusSeq{1,1} = seqconsensus(PoolDataStruc(i).ExampleSeqs(1:100,1));
        PoolDataStruc(i).ConsensusSeq{1,2} = seqconsensus(PoolDataStruc(i).ExampleSeqs(1:100,2));
        PoolDataStruc(i).LibFromAlignment = seqconsensus([PoolDataStruc(i).ConsensusSeq{1,1},seqrcomplement(PoolDataStruc(i).ConsensusSeq{1,2})]);
        disp(strcat(['Completed Analysis of Pool ID: ', num2str(i),'. Number of pools included: ',num2str(NumPools)]))
    end
    %{
    %Tells you what it found out
    for i=1:NumFiles
        tmpfilename = strcat(FastQDirectory, FastQFiles(i).name);
        for j=1:NumPools
            for k = 1:length(PoolDataStruc(j).FilePath)
                if strcmp(PoolDataStruc(j).FilePath(k),tmpfilename)
                    tmppoolname = PoolDataStruc(j).PoolName;
                    tmppoolnum = j;
                end
            end
        end
        disp(strcat(['File: ',num2str(i),': FastQName: ',FastQFiles(i).name, '  Pool Name: ', tmppoolname, 'PoolNumber: ', num2str(tmppoolnum)]))
    end
    %}
    
    
end