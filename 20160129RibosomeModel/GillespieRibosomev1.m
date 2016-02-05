clear all
close all
rand('state',sum(100*clock)); %#ok<RAND>

FullTrials = 10;

GeneArray = [5,10,15,20,25,30,35,40,45,50];
RunsArray = ones(50,1);%[20,20,20,20,20,20,20,20,20,20];%[120,60,40,30,24,20,17,14,13,12];
RunsArray = RunsArray .*30;

AllRibosomesPerGene = zeros(3,max(RunsArray),max(GeneArray));
for k = 1:FullTrials
    Case = k;
    disp(Case)

    

    NumGenes = GeneArray(Case);
    RibosomesInitial = 50*NumGenes;
    Runs = RunsArray(k);
    
    kB = .01;
    
    kON = .1;
    kOFF = 1;
    alpha = .01;
    
    tMax = 60;
    dt = 1;
    tspan = 0:dt:tMax;
    tspan2 = length(tspan);


    BurstDurationTrack = zeros(100000,Runs);
    BurstTimesTrack = zeros(100000,Runs);
    a = zeros(NumGenes*2,1);
    TraceVar = zeros(Runs,NumGenes);
    TraceAuto = zeros(tspan2,Runs*NumGenes);
    TraceAutoNorm = TraceAuto;

    Genes = zeros(1,NumGenes*2);
    Genes(NumGenes+1:end) = 1;
    BoundRibosomeArray = zeros(1000,1);
    mRNATotals = zeros(1,NumGenes);
    aGenes = zeros(NumGenes*2,1);
    
    GeneRuns = zeros(Runs,tspan2,NumGenes);
    RibosomesPerGene = zeros(Runs,NumGenes);
    
    %Build Rxn Matrix and Propensity for genes
    RxnMatrixGene = zeros(NumGenes*2,NumGenes*2);
    for i = 1:NumGenes
        RxnMatrixGene(i,i) = 1;
        RxnMatrixGene(i,NumGenes+i) = -1;
        aGenes(i) = kON;
    end
    for i = 1:NumGenes
        RxnMatrixGene(i+NumGenes,i) = -1;
        RxnMatrixGene(i+NumGenes,NumGenes+i) = 1;
        aGenes(i+NumGenes) = kOFF;
    end
    
    %
    
    for i = 1:Runs

        MaxOutput = 100000; %Maximum size expected of the output file
        %T = zeros(MaxOutput, 1); %Time tracking
        %X = zeros(1, NumSpecies); %Species number tracking
        %T(1)     = 0;
        %X(1,:)   = x0;
        RxnCount = 1;
        T = 0; %Time
        mRNAArray = zeros(1000,2);
        Ttrack = zeros(length(tspan), 1); %Time tracking
        XGenes = zeros(length(tspan), NumGenes*2); %Species number tracking
        XmRNA = zeros(length(tspan),length(mRNAArray(:,1)));
        XRibosomes = XmRNA;
        Ribosomes = RibosomesInitial;
        XGenes(1,:)   = Genes;
        xCurrentGenes = Genes;
        xCurrentmRNA = zeros(1000,1);
        xCurrentRibosomes = zeros(1000,1);
        mRNATotals = zeros(1,NumGenes);
        
        RecordTime = dt; %Recording time
        RecordCount = 2;
        OnTrack = 0;
        mRNAnum = 0;
        count = 1;
        OnDuration = zeros(MaxOutput,1);
        OnTimes = zeros(MaxOutput,1);
        %%%%%%%%%%%%%%%%%%%%%%
        %Gillespie Simulation%
        %%%%%%%%%%%%%%%%%%%%%%

        while T <= tMax

            % Calculate reaction propensities
            x = xCurrentGenes;

            for j = 1:NumGenes
                aGenes(j) = kON*xCurrentGenes(NumGenes+j);
                aGenes(j+NumGenes) = kOFF*xCurrentGenes(j);
                aAlpha(j) = xCurrentGenes(j)*alpha;
            end
            
            aRibosomes = xCurrentmRNA(:,1) .* Ribosomes .* kB;
            
            
            % Compute tau and mu using random variables
            a0 = sum(aGenes) + sum(aAlpha) + sum(aRibosomes);
            r = rand(1,2);
            %tau = -log(r(1))/a0;
            tau = (1/a0)*log(1/r(1));

            %Store information if at time
            tau0 = tau;
            if tau0 > dt
                while tau0 > dt
                    XGenes(RecordCount,:) = xCurrentGenes;
                    XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
                    XRibosomes(RecordCount,:) = xCurrentRibosomes;
                    RecordCount = RecordCount + 1;
                    RecordTime = RecordTime + dt;
                    tau0 = tau0 - dt;
                end
            end

            if T + tau > RecordTime
                XGenes(RecordCount,:) = xCurrentGenes;
                XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
                XRibosomes(RecordCount,:) = xCurrentRibosomes;
                RecordCount = RecordCount + 1;
                RecordTime = RecordTime + dt;

            end
            %[~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]);
            if r(2)*a0 < sum(aGenes)
                
                mu  = find((cumsum(aGenes) >= r(2)*a0),1,'first');
                T   = T  + tau;
                xCurrentGenes = xCurrentGenes + RxnMatrixGene(mu,:);
                
            elseif r(2)*a0 < sum(aGenes) + sum(aAlpha)
                mu  = find((sum(aGenes) + cumsum(aAlpha) >= r(2)*a0),1,'first');
                T   = T  + tau;
                mRNATotals(mu) = mRNATotals(mu) + 1;
                mRNAnum = mRNAnum + 1;
                xCurrentmRNA(mRNAnum,1) = 1;
                mRNAArray(mRNAnum,1) = 1;
                mRNAArray(mRNAnum,2) = mu;
                
            else
                mu  = find((sum(aGenes) + sum(aAlpha) + cumsum(aRibosomes) >= r(2)*a0),1,'first');
                %find the next change in state before tau
                T   = T  + tau;
                %xCurrentGenes = xCurrentGenes + RxnMatrixGene(mu,:);
                xCurrentRibosomes(mu) = xCurrentRibosomes(mu) + 1;
                Ribosomes = Ribosomes - 1;
            end
            
            

        end


        % Record output
         XGenes(RecordCount,:) = xCurrentGenes;
         XmRNA(RecordCount,:) = xCurrentmRNA(:,1);
         XRibosomes(RecordCount,:) = xCurrentRibosomes;
         
         mRNATotnum(i) = sum(mRNATotals);
         for j = 1:mRNATotnum(i)
            RibosomesPerGene(i,mRNAArray(j,2)) = RibosomesPerGene(i,mRNAArray(j,2)) + xCurrentRibosomes(j);
         end
         
        % Record output
        count = count - 1;

        GeneRuns(i,:,:) = XGenes(1:tspan2,1:NumGenes);
        mRNARuns(i,:,:) = XmRNA(1:tspan2,:);
        RibosomeRuns(i,:,:) = XRibosomes(1:tspan2,:);
        mRNAIdx(i,:) = mRNAArray(:,2); 



    end

   AllRibosomeRuns(k,:,:,:) = RibosomeRuns;
   AllRibosomesPerGene(k,1:Runs,1:NumGenes) = RibosomesPerGene;
   AllmRNATot(k,:) = mRNATotnum;
   AllmRNAIdx(k,:,:) = mRNAIdx;
   
end
save AllData
%%
%Calculations
for k = 1:FullTrials
    for j = 1:Runs
        temp = AllRibosomeRuns(k,j,:,:);
        temp = reshape(temp,[tspan2,1000]);
        tempsum = sum(temp,1);
        [RibosomeIdx, misc] = find(tempsum);
        for i = 1:length(RibosomeIdx)
            EndRibosomeCount(k,j,i) = temp(end,i);
            
        end
        tempsum = sum(temp,2);
        AllTotalRibosomes(k,j,:) = tempsum;
        EndmRNAwRibosomes(k,j) = length(RibosomeIdx); 
    end
    
    
    temp = mean(AllTotalRibosomes(k,:,:),2);
    TotalRibosomesGenTrend(k,:) = temp(:);
    
   
    for j = 1:Runs
        TotalRibosomeSS(k,j) = AllTotalRibosomes(k,j,end);
        temp = AllTotalRibosomes(k,j,:);
        temp2 = TotalRibosomesGenTrend(k,:);
        TotalRibosomeA(k,j,:) = temp(:) - temp2(:);
        TotalRibosomeVar(k,j) = var(TotalRibosomeA(k,j,:));
        TotalRibosomecv2(k,j) =  TotalRibosomeVar(k,j)/TotalRibosomeSS(k,j)^2;
        
    end
end
%%
%plot Total Ribosomes over time
c = colormap(jet(FullTrials));
hold on
for k = 1:FullTrials
    for j = 1:Runs
        temp = AllTotalRibosomes(k,j,:);
        plot(temp(:),'color',c(k,:))
    end
end
xlabel('Time','FontSize',15)
ylabel('Total Ribosomes','FontSize',15)
title('Total Ribosome over time')
saveas(gcf,'RibsomesTotalOverTime.jpg')
%%
%plot cv2 vs abundance
figure
hold on
c = colormap(jet(FullTrials));
for k = 1:FullTrials
    
    plot(TotalRibosomeSS(k,:),TotalRibosomecv2(k,:),'color',c(k,:),...
        'linestyle','none','marker','.','markersize',10)

end

for k = 1:FullTrials
    tempSS = TotalRibosomeSS(k,:);
    tempSS(tempSS == 0) = NaN;
    s = isnan(tempSS);
    SSmean = geomean(tempSS(~s));
    tempcv2 = TotalRibosomecv2(k,:);
    tempcv2(tempcv2 == Inf) = NaN;
    
    s = isnan(tempcv2);
    cv2mean = geomean(tempcv2(~s));
    plot(SSmean,cv2mean,'marker','o','markerfacecolor',c(k,:),...
        'linestyle','none','markersize',8,'markeredgecolor','k')

end

set(gca,'XScale','log');
set(gca,'YScale','log');
axis([.1 100000 10^-4 100]);
xlabel('Ribosome Abundance','FontSize',15)
ylabel('Ribosome cv^2','FontSize',15)
title('Total Ribosome cv^2 vs Abundance')
saveas(gcf,'Ribsomescv2abundance.jpg')
%%
figure
hold on
c = colormap(jet(FullTrials));
for k = 1:FullTrials
    for j = 1:Runs
        temp = AllRibosomeRuns(k,j,:,:);
        temp = reshape(temp,[tspan2,1000]);
        tempsum = sum(temp,1);
        [RibosomeIdx, misc] = find(tempsum);
        for i = 1:length(RibosomeIdx)
            EndRibosomeCount(k,j,i) = temp(end,i);
            
        end
        tempsum = sum(temp,2);
        AllTotalRibosomes(k,j,:) = tempsum;
        %[f,x] = ecdf(temp(end,1:length(RibosomeIdx)));
        %plot(x,f,'color',c(k,:))
        if isempty(RibosomeIdx)
        else
        plot(GeneArray(k),temp(end,1:length(RibosomeIdx)),'color',c(k,:),...
            'linestyle','none','marker','.','markersize',10)
        EndmRNAwRibosomes(k,j) = length(RibosomeIdx); 
        end
    end
end
xlabel('Concurrent Genes','FontSize',15)
ylabel('Ribosomes per Transcript','FontSize',15)
title('Number of Ribsomes per transcript')
saveas(gcf,'RibsomespermRNA.jpg')
%%
figure
hold on
for k = 1:FullTrials
    for j = 1:Runs
        temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
        temp = temp(:);
        plot(GeneArray(k),temp,'linestyle','none',...
            'marker','.','color',c(k,:),'markersize',10)
    end
end
xlabel('Concurrent Genes','FontSize',15)
ylabel('Ribosomes per Gene','FontSize',15)
title('Number of Ribsomes per Gene')
saveas(gcf,'RibsomesperGene.jpg')

%%
figure
hold on
c = colormap(jet(FullTrials));
for k = 1:FullTrials
    for j = 1:Runs
        temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
        temp = temp(:);
        temp(temp == 0) = [];
        meanBS(k,j) = mean(temp);
        meanBF(k,j) = length(temp);
        plot(meanBF(k,:), meanBS(k,:),'linestyle','none',...
            'marker','.','color',c(k,:),'markersize',12)
    end
end
for k = 1:FullTrials
    name = sprintf('%g Genes',GeneArray(k));
    tempBS = nanmean(meanBS(k,:));
    tempBF = nanmean(meanBF(k,:));
    linestore(k) = plot(tempBF, tempBS,'linestyle','none',...
        'marker','o','markerfacecolor',c(k,:),'markersize',8,...
        'markeredgecolor','k','displayname',name);

end
legend(linestore,'location','northeast')
% set(gca,'XScale','log');
% set(gca,'YScale','log');
ylabel('Mean Intensity ON Gene(BS)','FontSize',15)
xlabel('Mean Number ON (BF)','FontSize',15)
title('Burst Frequency Total Ribosome Per Gene ')
saveas(gcf,'BFvsBS.jpg')
%figure
%%
figure
hold on
c = colormap(hsv(Runs));
for k = 10
    for j = 1:Runs
        temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
        temp = temp(:);
        temp(temp == 0) = [];
        meanBS = mean(temp);
        meanBF = length(temp);
        plot(meanBF, meanBS,'linestyle','none',...
            'marker','.','color',c(j,:),'markersize',10)
    end
end
% set(gca,'XScale','log');
% set(gca,'YScale','log');
ylabel('Mean Intensity ON Gene(BS)','FontSize',15)
xlabel('Mean Number ON (BF)','FontSize',15)
title('Burst Frequency Total Ribosome Per Gene ')
%for a particular number of genes and run
c = colormap(hsv(Runs));
figure
hold on
for k = 10
    for j = 1:Runs
        %grab the total ribosomes per gene
        temp = AllRibosomesPerGene(k,j,1:GeneArray(k));
        %reshape to make pretty
        temp = temp(:);
        %find the corresponding mRNAs to this gene index
        temp2 = find(AllmRNAIdx(k,j,:)==14);
        temp3 = zeros(101,1);
        %for each mRNA
        for i = 1:length(temp2)
            %Find the time trace of the ribosomes to that mRNA
            temp4 = AllRibosomeRuns(k,j,:,temp2(i));
            %reshape
            temp4 = temp4(:);
            %sum them together to get a ribosome per gene trace
            temp3 = temp3 + temp4;
        end
        plot(temp3,'color',c(j,:))
    end
end


% hold on
% for i = 1:Runs
%     tempplot = AllRibosomeRuns(10,1,:,i);
%     plot(tempplot(:))
% end
% edges = 0:20:400;
% for k = 1:FullTrials
%     figure
%     temp = AllRibosomesPerGene(k,:,1:GeneArray(k));
%     temp = temp(:);
%     histogram(temp(:),edges,'normalization','probability')
%     xlabel('Ribosomes Per Gene','FontSize',15)
%     ylabel('Count','FontSize',15)
%     name = sprintf('Total Ribosomes for %g Genes',GeneArray(k));
%     title(name)
%     name = sprintf('HistogramTotRibosome%gGenes.jpg',GeneArray(k));
%     saveas(gcf,name)
% end
