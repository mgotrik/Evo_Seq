function [S,PoolDataStruct] = ES_GetFamilyFitness(PoolDataStruct,CurrThreshold,PrevThreshold)
%% Created 03/13/2019
%   Step 1: Loads a given round counted sequence information
%   Searches through the files under PoolDataStruct.ADS.FilePaths{4} Structure
%   with .Sequence and .count fields indicating the count of the sequence
%   in a given round.
%   Step 2: 
    
    
    FieldNamesInOrder = 
    
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
    
    %{
    Structure(SelectionNo).FamilyStruct(FamilyNo)
    
    
    %}                                                         
    
    ParentIDs = [];
    FamilyIDs
    SelectionCount = 0;
    for i = 1:NumFiles %Go through the files
        if exist(PoolDataStruct(i).ADS.FilePaths{4,1},'file') %This pool has been counted already
            if PoolDataStruct(i).ADS.Parent==0
                SelectionCount = SelectionCount+1;
                ParentIDs(SelectionCount) = i;
                FamilyIDs(SelectionCount) = PoolDataStruct(i).ADS.Family;
            end
        end
    end
            
    
    %SS(SelectionID).PoolStruct(PoolID)
    for SelectionID = 1:SelectionCount %Going through each selection
        Parent = ParentIDs(SelectionID);
        Family = FamilyIDs(SelectionCount);
        PoolStruct = struct();
        if length(Family)<2 %
            continue
        end
        for PoolID = fliplr(Family(2:end)) %Going through each pool
    %
            tic
            disp(strcat("Beginning in file", num2str(PoolID)))  %
            PoolDataStruct(PoolID).ADS.FilePaths{5,1} = strrep(PoolDataStruct(i).ADS.FilePaths{4},"_counted","_fitness"); %
            CurrFieldName = PoolDataStruct(PoolID).PoolName{1};

            if exist(PoolDataStruct(PoolID).ADS.FilePaths{5},"file")==2
                S.(CurrFieldName)=load(PoolDataStruct(PoolID).ADS.FilePaths{5,1},CurrFieldName);
                S.(CurrFieldName)=S.(CurrFieldName).(CurrFieldName);
            else
                % MAY NEED TO REFACTOR THIS WITH REAL VALUES
                S.(CurrFieldName)=load(PoolDataStruct(PoolID).ADS.FilePaths{4,1},'Sequences');
                S.(CurrFieldName) = S.(PoolDataStruct(PoolID).PoolName{1}).Sequences;
            end

            if ~isfield(S.(CurrFieldName),'fitness')
                S.(CurrFieldName)(1).fitness = 1;
            end
            for PrevIndex = fliplr(Family(1):find(Family==PoolID)-1)
                PoolDataStruct(PrevIndex).ADS.FilePaths{5,1} = strrep(PoolDataStruct(PrevIndex).ADS.FilePaths{4,1},"_counted","_fitness");
                PrevFieldName = PoolDataStruct(PrevIndex).PoolName{1};%fieldnames(fieldnames(S.(PoolDataStruct(tmp).PoolName{1}))){1}
                % WILL NEED TO REFACTOR THIS TO MAKE IT WORK
                if exist(PoolDataStruct(PrevIndex).ADS.FilePaths{5},"file")==2
                    S.(PrevFieldName)=load(PoolDataStruct(PrevIndex).ADS.FilePaths{5},'Sequences');
                    S.(PrevFieldName)=S.(PrevFieldName).Sequences;
                else
                    S.(PrevFieldName)=load(PoolDataStruct(PrevIndex).ADS.FilePaths{4});
                    S.(PrevFieldName) = S.(PrevFieldName).(PoolDataStruct(PrevIndex).ADS.FilePaths{4}StructName);
                end
                if ~isfield(S.(PrevFieldName),'fitness')
                    S.(PrevFieldName)(1).fitness = 1;
                end
                [S.(PrevFieldName),S.(CurrFieldName)] = ComparePoolsForFitness(S.(PrevFieldName),S.(CurrFieldName),PrevIndex,PoolID,PoolDataStruct,CurrThreshold,PrevThreshold);
                tempStruct = struct;
                if PrevIndex==1 && PoolID==2 %If this is the last one we are looking at
                    S.(PrevFieldName) = orderfields(S.(PrevFieldName),FieldNamesInOrder);
                end
                tempStruct.(PrevFieldName) = S.(PrevFieldName);
                S = rmfield(S,PrevFieldName);
                save(PoolDataStruct(PrevIndex).ADS.FilePaths{5},'-struct', 'tempStruct');
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
            save(PoolDataStruct(PoolID).ADS.FilePaths{5},'-struct', 'tempStruct');
            clearvars tempStruct
        end
    end

end

    %%
