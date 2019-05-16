function [PoolDataStruct] = ES_getPDSLibInfo(PoolDataStruct,ParentPools,LibIDs)

[LibStruct] = excel2struct('Libraries.xls',[]);
struct2table(LibStruct)

for i = 1:length(LibStruct)
    LibStruct(i).Sequence = upper(LibStruct(i).Sequence)
    LibStruct(i).FP = upper(LibStruct(i).FP);
    LibStruct(i).RP = upper(LibStruct(i).RP);
end

NumParents = 0; %Counting variables for the total number of different pools/selections happening

for i = 1:length(ParentPools) %Iterates through each provided [pool,parent] element
    CurrPool =  ParentPools(i,1);
    ParPool = ParentPools(i,2);
    if ParPool==0 %|| i==length(ParentPools) %If this is a new pool family/selection
        NumParents = NumParents+1; 
        FamilyID{NumParents} = CurrPool; %Assigns it a new family ID
    else %If this is not a new pool
        for j = 1:length(FamilyID) %It checks all Families
            if find(FamilyID{j}==ParPool) %To see if its parent is in a family, and updates family info
                FamilyID{j} = [FamilyID{j}, CurrPool];
                break
            end
        end
    end
end

PoolCount = 0; %Counting variables for the number of unique selections/pools

for i = 1:length(FamilyID) %Iterates through each family
    for j = FamilyID{i} %Iterates through each member
        PoolCount = PoolCount+1; %The number of pools we have seen so far
        PoolDataStruct(PoolCount).Parent = ParentPools(PoolCount,2);
        PoolDataStruct(PoolCount).Family = FamilyID{i};
        PoolDataStruct(j).LibNames = {LibStruct(LibIDs{i}).Name}';
        PoolDataStruct(j).LibStruct = LibStruct(LibIDs{i});
    end
end

%% Works up until here:

for i = 1:PoolCount
    if isempty(PoolDataStruct(i).ReadLen)
        continue
    end
    Libraries = PoolDataStruct(i).LibStruct;

    
    for j = 1:length(Libraries) %Iterates through each library in the read
        % Begins with the full inserted sequence
        TempLib = Libraries(j).Sequence;
        % Splits where the TrueSeq Read1 primer is (or whatever primer
        % where R1 sequencing begins
        TempLib = strsplit(TempLib,Libraries(j).Read1);
        
        % Splits where the TrueSeq Read2 primer is (or whatever primer
        % where R2 sequencing begins
        TempLib = strsplit(TempLib{2},seqrcomplement(Libraries(j).Read2));
        % Assigns this as the sequenced library
        SeqLib = TempLib{1};
        LibLength = length(SeqLib);
        Libraries(j).SeqLib= SeqLib;
        
        % Finds the maximum length R1 can be (if sequenced library <= total
        % R1 length, index 2 is the minimum of these two)
        R1MaxLen = min(PoolDataStruct(i).ReadLen(1),length(SeqLib));
        R1Index = [1, R1MaxLen];
        SeqLibR1Idx = [R1Index(1):R1Index(2)];
        
        % Finds the maximum length R2 can be (if sequenced library <= total
        % R2 length, index 1 is the maximum of these two)
        R2MinLen = max(1,LibLength - PoolDataStruct(i).ReadLen(2)+1);
        R2Index = [R2MinLen, LibLength];
        SeqLibR2Idx = [R2Index(1):R2Index(2)];
        
        %Identifies the overlap region as the region between R2Min and
        %R1Max
        OLIdx = [R2MinLen, R1MaxLen];
        SeqLibOLIdx = [R2MinLen:R1MaxLen];
        
        Libraries(j).R1Index = R1Index;
        Libraries(j).R2Index = R2Index;
        Libraries(j).OLIdx = OLIdx;
        Libraries(j).LibRead1 = SeqLib(SeqLibR1Idx);
        Libraries(j).LibRead2 = SeqLib(SeqLibR2Idx);
        Libraries(j).LibRead2rc = seqrcomplement(Libraries(j).LibRead2);
        Libraries(j).ReadOL = SeqLib(SeqLibOLIdx);
        
        FP = Libraries(j).FP;
        RPc = seqrcomplement(Libraries(j).RP);
        
        
        %{
        if isempty(strfind(PDS(i).Seq_Read1{j},ForwardPrimer))
            ReadFlag = 2; %Consider only Read 2
        elseif isempty(strfind(PDS(i).Seq_Read2{j},RPc))
            ReadFlag = 1; %Consider only Read 1
        elseif ~isempty(strfind(PDS(i).Seq_Read1{j},ForwardPrimer)) && ~isempty(strfind(PDS(i).Seq_Read2{j},RPc))
        %}
        
        if ~isempty(strfind(Libraries(j).LibRead1,RPc))
            Libraries(j).ReadFlag = 1; %Consider only Read 1
        elseif ~isempty(strfind(Libraries(j).LibRead2,FP))
            Libraries(j).ReadFlag = 2; %Consider only Read 2
        elseif ~isempty(strfind(Libraries(j).LibRead1,FP)) && ~isempty(strfind(Libraries(j).LibRead2,RPc))
            Libraries(j).ReadFlag = 3; %Library is covered in both Read1 and Read2
        else
            error('Logic failed, oh no!')
        end
        
    end
    PoolDataStruct(i).LibStruct = Libraries;    
end
%}

end