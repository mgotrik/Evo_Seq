function [S,PoolDataStruct] = ArrangePoolsForSequenceAnalysis(PoolDataStruct,CurrThreshold,PrevThreshold,FileOrder)
%% Created 03/13/2019
%   Step 1: Loads a given round counted sequence information
%   Searches through the files under PoolDataStruct.SavedMatFile Structure
%   with .Sequence and .count fields indicating the count of the sequence
%   in a given round.
%   Step 2: 


    PoolNoStart = 4;

    NumFiles = size(PoolDataStruct,2);
    %FileNames(1:7) = {char(strcat(PoolDataStruct(1:7).FileDir,PoolDataStruct(1:7).PoolName,'_counted.mat'))}
    %S.(FileNames{1:2}) = load(FileNames(1:2))
    S = struct;
    
    %S.(PoolDataStruct(NumFiles).PoolName{1})=load(AllFileNames{NumFiles});
    %CurrRoundSequences = S.(PoolDataStruct(NumFiles).PoolName)
    %FileNames = {'Seqs_BFAR1','Seqs_BFAR2','Seqs_BFAR3','Seqs_BFAR3ep','Seqs_BFAR4',... 
    %    'Seqs_BFAR4H','Seqs_BFAR4L'};
    FieldNamesRCN = cell(0,0);
    FieldNamesIdx = cell(0,0);
    for i=PoolNoStart:NumFiles
        CurrRPoolName = PoolDataStruct(i).PoolName;
        FieldNameCurrRCN = cellstr(strcat(CurrRPoolName,"_RPM"));
        FieldNameCurrRIdx = cellstr(strcat(CurrRPoolName,"_Idx"));
        FieldNamesRCN = [FieldNamesRCN,FieldNameCurrRCN];
        FieldNamesIdx = [FieldNamesIdx,FieldNameCurrRIdx];
    end
    FieldNamesInOrder = [{'Sequence'},{'count'},{'fitness'},FieldNamesIdx,FieldNamesRCN];
    disp(FieldNamesInOrder)
    for CurrIndex=fliplr(PoolNoStart+1:NumFiles)
        %
        tic
        disp(strcat("Beginning in file", num2str(CurrIndex)))
        PoolDataStruct(CurrIndex).SavedMatFileFitness = strrep(PoolDataStruct(CurrIndex).SavedMatFile,"_counted","_fitness");
        CurrFieldName = PoolDataStruct(CurrIndex).PoolName{1};
        
        if exist(PoolDataStruct(CurrIndex).SavedMatFileFitness,"file")==2
            S.(CurrFieldName)=load(PoolDataStruct(CurrIndex).SavedMatFileFitness,CurrFieldName);
            S.(CurrFieldName)=S.(CurrFieldName).(CurrFieldName);
        else
            S.(CurrFieldName)=load(PoolDataStruct(CurrIndex).SavedMatFile,PoolDataStruct(CurrIndex).SavedMatFileStructName);
            S.(CurrFieldName) = S.(PoolDataStruct(CurrIndex).PoolName{1}).(PoolDataStruct(CurrIndex).SavedMatFileStructName);
        end
        
        if ~isfield(S.(CurrFieldName),'fitness')
            S.(CurrFieldName)(1).fitness = 1;
        end
        for PrevIndex = fliplr(PoolNoStart:CurrIndex-1)
            PoolDataStruct(PrevIndex).SavedMatFileFitness = strrep(PoolDataStruct(PrevIndex).SavedMatFile,"_counted","_fitness");
            PrevFieldName = PoolDataStruct(PrevIndex).PoolName{1};%fieldnames(fieldnames(S.(PoolDataStruct(tmp).PoolName{1}))){1}
            
            if exist(PoolDataStruct(PrevIndex).SavedMatFileFitness,"file")==2
                S.(PrevFieldName)=load(PoolDataStruct(PrevIndex).SavedMatFileFitness,PrevFieldName);
                S.(PrevFieldName)=S.(PrevFieldName).(PrevFieldName);
            else
                S.(PrevFieldName)=load(PoolDataStruct(PrevIndex).SavedMatFile);
                S.(PrevFieldName) = S.(PrevFieldName).(PoolDataStruct(PrevIndex).SavedMatFileStructName);
            end
            if ~isfield(S.(PrevFieldName),'fitness')
                S.(PrevFieldName)(1).fitness = 1;
            end
            [S.(PrevFieldName),S.(CurrFieldName)] = ComparePoolsForFitness(S.(PrevFieldName),S.(CurrFieldName),PrevIndex,CurrIndex,PoolDataStruct,CurrThreshold,PrevThreshold);
            tempStruct = struct;
            if PrevIndex==1 && CurrIndex==2 %If this is the last one we are looking at
                S.(PrevFieldName) = orderfields(S.(PrevFieldName),FieldNamesInOrder);
            end
            tempStruct.(PrevFieldName) = S.(PrevFieldName);
            S = rmfield(S,PrevFieldName);
            save(PoolDataStruct(PrevIndex).SavedMatFileFitness,'-struct', 'tempStruct');
            clearvars tempStruct
            disp(strcat("Finished comparing ",CurrFieldName," to ",PrevFieldName))
        end
        toc
        disp(strcat("Finished comparing ",CurrFieldName," to all other files"))
        tempStruct = struct;
        %Sorts S such that it goes Indx, by round, RPM, by round, in order:
        S.(CurrFieldName) = orderfields(S.(CurrFieldName),FieldNamesInOrder);
        
        tempStruct.(CurrFieldName) = S.(CurrFieldName);
        S = rmfield(S,CurrFieldName);
        save(PoolDataStruct(CurrIndex).SavedMatFileFitness,'-struct', 'tempStruct');
        clearvars tempStruct
    end
    
end

    %%
