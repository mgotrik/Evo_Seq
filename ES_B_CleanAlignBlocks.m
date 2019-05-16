
function [CombSeqs,AlignFail,GappedReads] = ES_B_CleanAlignBlocks(TempSeqF,TempSeqR,BlockSize)
    %%
    AlignFail = 0;
    GappedReads = 0;
    
    NewSeq = cell(BlockSize,1);
    NewQual = cell(BlockSize,1);
    NewHead = cell(BlockSize,1);
    for k = 1:BlockSize %Iterates through each sequence in the block
        
        %Gets the sequence and quality
        SeqF = TempSeqF(k).Sequence; 
        SeqR = TempSeqR(k).Sequence;
        QualF = TempSeqF(k).Quality;
        QualR = TempSeqR(k).Quality;
        
        

        %Aligns the two with an extensive penalty for a gap (gaps shouldn't
        %exist that often when aligning the reads, but sometimes do)
        
        %29.6s w 2212 fail; 29 gapped out of 30k:
        %[~,b,c]=swalign(SeqF,SeqR,'GapOpen',15,'ExtendGap',10,'Alphabet','NT');
        
        %28.9s w 2212 fail; 29 gapped out of 30k:
        [~,b,c]=swalign(SeqF,SeqR,'GapOpen',15,'Alphabet','NT');
        
        LenOL = size(b,2);
        %32.8s w 2212 fail; 29 gapped out of 30k:
        %{ 
        align=localalign(SeqF,SeqR,'GapOpen',15,'Alphabet','NT');
        b = align.Alignment{1};
        c = align.Start;
        %}
        
        if  contains(b(2,:),'||||||||') %LenOL>5 %contains(b(2,:),'|||||||') %find(b(2,:)=='|||||||') %LenOL>5 %arbitrary >5

            %LenF = LenOL - count(b(1,:),'-'); %How many bases in OL of F
            %LenR = LenOL - count(b(3,:),'-'); %How many bases in OL of R
            
            if (count([b(1,:),b(3,:)],'-'))
                GappedReads = GappedReads+1;
                AlignFail = AlignFail+1;
                %{
                if sum(QualBR)/length(QualBR) > sum(QualBF)/length(QualBF)
                    SeqB = SeqBR;
                    QualB = QualBR;
                else
                    SeqB = SeqBR;
                    QualB = QualBR;
                end
                %}
                %GappedSeq(GappedReads) = {[SeqA,SeqB,SeqC]};
                %GappedQual(GappedReads) = {[QualA,QualB,QualC]};
            else
                
                SeqA = SeqF(1:c(1)-1);
                SeqC = SeqR(c(2)+size(b,2):end);
                QualA = QualF(1:c(1)-1);
                QualC = QualR(c(2)+size(b,2):end);

                SeqBF = SeqF(c(1):c(1)+LenOL-1);
                SeqBR = SeqR(c(2):c(2)+LenOL-1);
                QualBF = QualF(c(1):c(1)+LenOL-1);
                QualBR = QualR(c(2):c(2)+LenOL-1);
                
                QualMat = [QualBF;QualBR];
                SeqMat = [SeqBF;SeqBR];

                [~,SortIdx] = sort(QualMat,1,'descend');
                SeqMat = sortrows(SeqMat,SortIdx(1,:));
                QualMat = sortrows(QualMat,SortIdx(1,:));
                for i = 1:length(SortIdx)
                    SeqB(i) = SeqMat(SortIdx(1,i),i);
                    QualB(i) = QualMat(SortIdx(1,i),i);
                end
                NewHead(k-AlignFail) = {TempSeqF(k).Header};
                NewSeq(k-AlignFail) = {[SeqA,SeqB,SeqC]};
                NewQual(k-AlignFail) = {[QualA,QualB,QualC]};
            end

            clearvars SeqB QualB
        else
            AlignFail = AlignFail+1;
        end
        
    end
    NewHead(BlockSize-AlignFail+1:end)='';
    NewSeq(BlockSize-AlignFail+1:end)='';
    NewQual(BlockSize-AlignFail+1:end)='';
    CombSeqs = struct('Header',NewHead,'Sequence',NewSeq,'Quality',NewQual);
    %CombSeqs = struct('Header',NewHead{:},'Sequence',NewSeq{:},'Quality',NewQual{:});
    %[CombSeqs.Header] = NewHead{:};
    %[CombSeqs.Sequence] = NewSeq{:};
    %[CombSeqs.Quality] = NewQual{:}     
    
end