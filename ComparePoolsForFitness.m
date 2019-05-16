function [EarlyRoundSeqs,CurrRoundSeqs] = ComparePoolsForFitness(EarlyRoundSeqs,CurrRoundSeqs, PrevPDSIndex,CurrPDSIndex, PoolDataStruct, CurrThreshold,PrevThreshold)
%{
    ParentSeqs = S.BFAR3ep
    CurrSeqs = S.BFAR4L
    ParentPDSIndex = 4
    CurrPDSIndex = 7
    CurrThreshold = 100
    ParentThreshold = 10
%}
    if PoolDataStruct(CurrPDSIndex).Parent==PrevPDSIndex
        isParentFlag = 1;
    else
        isParentFlag = 0;
    end


    PrevRCopyNo = PoolDataStruct(PrevPDSIndex).FilteredEntries;
    CurrRCopyNo = PoolDataStruct(CurrPDSIndex).FilteredEntries;
    PrevRNormalFactor = 1000000/PrevRCopyNo;
    CurrRNormalFactor = 1000000/CurrRCopyNo;
    
    
    CurrRPoolName = PoolDataStruct(CurrPDSIndex).PoolName;
    PrevRPoolName = PoolDataStruct(PrevPDSIndex).PoolName;
    FieldNameCurrRCN = strcat(CurrRPoolName,"_RPM");
    FieldNameCurrRIdx = strcat(CurrRPoolName,"_Idx");
    FieldNamePrevRCN =  strcat(PrevRPoolName,"_RPM");
    FieldNamePrevRIdx =   strcat(PrevRPoolName,"_Idx");
    %containers.Map(keyset,valueset)
    %keys = the individual sequences
    %values = the returned values; #1-CurrRindex; #2-CurrRCount; #3-
    %SEQMAP = containers.Map(CurrRSequences(1).Sequence,[1,CurrRSequences(1).count, 0, 0, 0]);
    
    SEQMAP = containers.Map;
    
    for i = 1:size(CurrRoundSeqs,1)
        SEQMAP(CurrRoundSeqs(i).Sequence) = [i,CurrRoundSeqs(i).count];
        if CurrRoundSeqs(i).count <CurrThreshold
            break
        end
    end
    disp(strcat('Finished mapping:',{' '},CurrRPoolName))
    for i = 1:size(EarlyRoundSeqs,1)
        TempSeq = EarlyRoundSeqs(i).Sequence;
        TempCount = EarlyRoundSeqs(i).count;
        if TempCount<PrevThreshold
            break
        end
        if isKey(SEQMAP,TempSeq)
            CurrValues = SEQMAP(TempSeq);
            CurrRIndex = CurrValues(1);
            CurrRCount = CurrValues(2);
            PrevRIndex = i;
            PrevRCount = TempCount;
            CurrRPM = round(CurrRCount*CurrRNormalFactor,2);
            PrevRPM = round(PrevRCount*PrevRNormalFactor,2);
            
            if isParentFlag==1
                CurrRoundSeqs(CurrRIndex).(FieldNameCurrRIdx) = CurrRIndex;
                CurrRoundSeqs(CurrRIndex).(FieldNameCurrRCN) = CurrRPM;
                CurrRoundSeqs(CurrRIndex).fitness = CurrRCount/PrevRCount;
                EarlyRoundSeqs(PrevRIndex).(FieldNamePrevRIdx) = PrevRIndex;
                EarlyRoundSeqs(PrevRIndex).(FieldNamePrevRCN) = PrevRPM;
            end
            CurrRoundSeqs(CurrRIndex).(FieldNamePrevRIdx) = PrevRIndex;
            CurrRoundSeqs(CurrRIndex).(FieldNamePrevRCN) = PrevRPM;
            EarlyRoundSeqs(PrevRIndex).(FieldNameCurrRIdx) = CurrRIndex;
            EarlyRoundSeqs(PrevRIndex).(FieldNameCurrRCN) = CurrRPM;
        end
    end
    
    disp(strcat('Finished comparing:',{' '},PrevRPoolName))
end