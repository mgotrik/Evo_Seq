function PoolDataStruct = ES_GetAnalysisInfo(PoolDataStruct,BlockSize,MinQuality,MaxPercentLowQualityBases,FilterDir)
%%
for PoolID = 1:length(PoolDataStruct)
    if isempty(PoolDataStruct(PoolID).FileDir)
        continue
    end

    if ~isempty(strfind(PoolDataStruct(PoolID).FileDir,'\'))
        SepChar = '\';
    else
        SepChar = '/';
    end
        

        
    for k = 1:2 %Iterating through the two read files
        StartPath = [PoolDataStruct(PoolID).FastQInfo(k).FilePath,SepChar];
        StartName = PoolDataStruct(PoolID).FastQInfo(k).Filename;
        FilterPath = [StartPath, 'FilteredFastQ',SepChar];
        FilterName = strrep(StartName,'.fastq','_filtered.fastq');
        FileNames{k} = [StartPath, StartName];
        FilterNames{k} = [FilterPath, FilterName];
    end
    Descriptor = {{'Starting FastQ files';'FastQ After seqfilter'}};
    
    %{
    'FastQNames',{PoolDataStruct(i).FastQInfo.Filename}
    'FilterDir',{FilterDir},'FilePaths',[FileNames;FilterNames],...
        'FilePathDescription',Descriptor,'FailedQualityCheck',0,...
        'MinQuality',MinQuality,'MaxPctLowQual',MaxPercentLowQualityBases)
    %}
    
    ADS = struct('Date',date,'BlockSize',BlockSize,'NumBlocks',0,'LastBlockSize',0,...
        'FilterDir',{FilterDir},'FilePaths',{[FileNames;FilterNames]},...
        'FilePathDescription',Descriptor,'FailedQualityCheck',0,...
        'MinQuality',MinQuality,'MaxPctLowQual',MaxPercentLowQualityBases,...
        'FilteredEntries',0);

    FilterFileName1 = ADS.FilePaths{2,1};    
    FilterFileName2 = ADS.FilePaths{2,2};
    
    if ~exist(ADS.FilterDir) %Makes directory if it doesn't exist
        mkdir(ADS.FilterDir);
    end
    PoolDataStruct(PoolID).ADS = ADS;
end
%%

tic
disp('Beginnning')


parfor i=1:PoolID
    if isempty(PoolDataStruct(i).FileDir)
        continue
    end
    if ~exist(PoolDataStruct(i).ADS.FilePaths{2,1},'file') || ~exist(PoolDataStruct(i).ADS.FilePaths{2,2},'file')
    PoolDataStruct(i).ADS.FilePathFilter = '';
    [PoolDataStruct(i).ADS.FilePathFilter,PoolDataStruct(i).ADS.FilteredEntries,PoolDataStruct(i).ADS.FailedQualityCheck] = seqfilter([{char(PoolDataStruct(i).ADS.FilePaths(1,1))},{char(PoolDataStruct(i).ADS.FilePaths(1,2))}],'PairedFiles',true,...
                             'Method','MaxPercentLowQualityBases',...
                             'Threshold',[MinQuality MaxPercentLowQualityBases],...
                             'OutputDir', FilterDir,...
                             'UseParallel', true,...
                             'WriteSingleton',false); %Does not write ones that have one read fail but one good
                         
    PoolDataStruct(i).ADS.FilteredEntries = PoolDataStruct(i).ADS.FilteredEntries(1);
    PoolDataStruct(i).ADS.FailedQualityCheck=PoolDataStruct(i).ADS.FailedQualityCheck(1);
    else
        InfoStruct = fastqinfo(PoolDataStruct(i).ADS.FilePaths{2,1})
        PoolDataStruct(i).ADS.FilteredEntries = InfoStruct.NumberOfEntries;
        PoolDataStruct(i).ADS.FailedQualityCheck = PoolDataStruct(i).FastQInfo(1).NumberOfEntries - InfoStruct.NumberOfEntries;
    end
    PoolDataStruct(i).ADS.NoSeqBlocks = ceil(PoolDataStruct(i).ADS.FilteredEntries(1)/PoolDataStruct(i).ADS.BlockSize(1)); %NoSeqBlocks = floor(SeqNoStart(i)/SeqBlockSize);
    PoolDataStruct(i).ADS.LastBlockSize = mod(PoolDataStruct(i).ADS.FilteredEntries(1),PoolDataStruct(i).ADS.BlockSize(1));
    NumFailed = PoolDataStruct(i).ADS.FailedQualityCheck(1)/PoolDataStruct(i).ADS.FilteredEntries(1);
    disp(['Number of low quality reads removed (%): ',num2str(PoolDataStruct(i).ADS.FailedQualityCheck(1)),' (',num2str(NumFailed*100),'). Pool:',PoolDataStruct(i).PoolName{1}]);
end
toc


end
