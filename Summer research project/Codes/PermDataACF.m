%Step 1: Obtain matrix of permeabilities for odd layers, with each row 
%representing one realisation or gap.

boxsize=1; %This is a box size of 5 pixels
fnames = dir('By column/P1L*_col*.txt'); %Directory for subsets of odd layers
totalboxes=floor((506-30)/boxsize); %Defines the number of whole boxes on interval
perm=zeros(408,476); %Initialise perms as a matrix with 17(gaps)x3(subsets)x8(layers) rows and 476 (totalboxes) columns
LayerNumber=1; %To assign permeability to correct row later
for j=2:7:51 %Each j represents the first valid subset (column 2) of an odd layer
    RowSpec=51*(LayerNumber-1);
    SubsetNumber=1;
    for s=1:3
        if s==1
           panelread=dlmread(fnames(j).name,' ',3,0); %Read first subset of layer from 4th row
        elseif s==2
           panelread=dlmread(fnames(j+2).name,' ',3,0); %Read third subset of layer from 4th row
        else 
           panelread=dlmread(fnames(j+4).name,' ',3,0); %Read fifth subset of layer from 4th row
        end
        for i=1:17 %i represents the number of the gap, we have at least 17 gaps in each layer
            startpos=30; %Reset the starting distance (this is 1480 pixels)
            for k=1:476 %k represents the box number along the gap
                panel=panelread; %Value reset, dlmread not called every time
                panel(:,size(panel,2))=[]; %Delete final column
                position=[startpos startpos+1]; %Vector specifying the start row and end row for each box
                panel=panel(position(1):position(2),:); %Take data in box
                panel=panel*35.3e-6; %Conversion to meters
                index=17+16*(i-1); %Find start of gap
                gapstart=panel(:,index); %Gives positions of gap start
                gapend=panel(:,index+1); %Gives positions of gap end
                gapsize=gapend-gapstart; %Vector of gap sizes
                gap_avg=mean(gapsize); %Find average gap size
        
                if gap_avg>0
                   width = 6.38e-3*8+gap_avg;
                   porosity = 1-(0.7*(6.35e-3*8)/(width));
                   if boxsize==0
                      orientation=0;
                   else
                      midstart=(gapstart(1)+gapend(1))/2;
                      midend=(gapstart(length(gapstart))+gapend(length(gapend)))/2;
                      orientation = atan(abs(midstart-midend)/(boxsize*35/1e-6))*180/pi;
                   end
                   [K1, K12, K2, K3] = calc_permeability_v2(1,porosity,orientation, 0,gap_avg,width,0.16e-3,3.5e-6);
                   perm(i+RowSpec+(SubsetNumber-1)*17,k)=K1; %Assign the permeability to the appropriate element in perm
                end
                startpos=startpos+1; %Move the focus along to the next box
            end
        end
        SubsetNumber=SubsetNumber+1;
    end
    LayerNumber=LayerNumber+1; %Add 1 here as a layer has been completed
end

%Step 2: Repeat above for the even layers (except layers 12 and 14)
enames=dir('By column/P1L*_row*.txt'); %Directory for subsets of even layers
permEven=zeros(408,476); %17(gaps per subset)*4(subsets per layer)*6(even layers)=408 rows 
LayerNumber=1;
for j=2:12:86
    if j==62||j==74 %Skip layers 12 and 14
        continue
    else
        RowSpec=68*(LayerNumber-1);
        SubsetNumber=1;
        for s=1:4 %Number of subsets we're interested in
            if s==1
                panelread=dlmread(enames(j).name,' ',3,0); %First subset
            elseif s==2
                panelread=dlmread(enames(j+3).name,' ',3,0); %Fourth subset
            elseif s==3
                panelread=dlmread(enames(j+6).name,' ',3,0); %Seventh subset
            else 
                panelread=dlmread(enames(j+9).name,' ',3,0); %Tenth subset
            end
            for i=1:17 %i represents the number of the gap, we have at least 17 gaps in each layer
                startpos=30; %Reset the starting distance (this is 1480 pixels)
                for k=1:476 %k represents the box number along the gap
                    panel=panelread; %Value reset, dlmread not called every time
                    panel(:,size(panel,2))=[]; %Delete final column
                    position=[startpos startpos+1]; %Vector specifying the start row and end row for each box
                    panel=panel(position(1):position(2),:); %Take data in box
                    panel=panel*35.3e-6; %Conversion to meters
                    index=17+16*(i-1); %Find start of gap
                    gapstart=panel(:,index); %Gives positions of gap start
                    gapend=panel(:,index+1); %Gives positions of gap end
                    gapsize=gapend-gapstart; %Vector of gap sizes
                    gap_avg=mean(gapsize); %Find average gap size
        
                    if gap_avg>0
                       width = 6.38e-3*8+gap_avg;
                       porosity = 1-(0.7*(6.35e-3*8)/(width));
                       if boxsize==0
                          orientation=0;
                       else
                          midstart=(gapstart(1)+gapend(1))/2;
                          midend=(gapstart(length(gapstart))+gapend(length(gapend)))/2;
                          orientation = atan(abs(midstart-midend)/(boxsize*35/1e-6))*180/pi;
                       end
                       [K1, K12, K2, K3] = calc_permeability_v2(1,porosity,orientation, 0,gap_avg,width,0.16e-3,3.5e-6);
                       permEven(i+RowSpec+(SubsetNumber-1)*17,k)=K1; %Assign the permeability to the appropriate element in perm
                    end
                    startpos=startpos+1; %Move the focus along to the next box
                end
            end
        SubsetNumber=SubsetNumber+1;
        end 
    LayerNumber=LayerNumber+1; %Add 1 here as a layer has been completed
    end
end

%Step 3: Combine the above two permeability matrices into one (408+408)x476
%matrix to be used later.

permeability=[perm;permEven];
%mu=mean(permeability,'all'); %TRIAL for combined mean - this completely
%changes shape of resulting acf and power spectral
%muVEC=mean(permeability,2); %Vector of means of each row (or realisation)

%Step 4: Find sample autocorrelation for each gap and plot along with
%average and confidence interval

lag=0:1:400; %Specify a desired vector of lags
acf=zeros(816,length(lag)); %Initialise matrix of autocorrelations
PSD=zeros(816,length(lag));
n=length(lag);
for k=1:816
    mu=mean(permeability(k,:));
    [acf1,~]=sampleacf(permeability(k,:),400,mu); %Call sampleacf.m to find sample autocorrelation over each gap
    acf(k,:)=acf1;
    PSD(k,:)=abs(fftshift(fft(acf1)));
end
lower=tinv(0.025,815); %Next two lines define the critical values of the appropriate t-distribution
upper=tinv(0.975,815);
Sn=sum(acf,1); %Sum columns of both matrices to get sum of independent realisations
SnPSD=sum(PSD,1);
ACFavg=mean(acf); %Find mean of each column in acf
PSDavg=mean(PSD);
sd=std(acf); %Find sample standard deviation of each column in acf
sdPSD=std(PSD); %Do same for PSD
lowerCI=Sn/816+(lower/sqrt(816))*sd; %Next lines define the lower and upper limits of the confidence intervals
upperCI=Sn/816+(upper/sqrt(816))*sd;
lowerCIPSD=SnPSD/816+(lower/sqrt(816))*sdPSD;
upperCIPSD=SnPSD/816+(upper/sqrt(816))*sdPSD;

figure(1);
dist=lag*5*0.0353; %Convert lag to mm
hold on %Rest of this code produces the plot
for l=1:10
    plot(dist,acf(l,:),'-','LineWidth',0.1);
end
plot(dist,ACFavg,'k-','LineWidth',2.5);
plot(dist,lowerCI,'--r','LineWidth',2.5);
plot(dist,upperCI,'--r','LineWidth',2.5);
hold off
ylim([-0.3 1.05]);
xlim([0 70]);
xlabel('Distance (mm)');
ylabel('Autocorrelation');
title('Sample autocorrelation (data) with 95% interval');

%% Section out 
figure(2);
subplot(2,1,1);
hold on
fs=10; %Sampling frequency
Ts=1/fs; %Sampling interval
df=(1/n)/Ts;
f=(-n/2:n/2-1)*df; %Define frequency domain
f=f-f(floor(n/2)+1)*ones(1,length(f)); %Centre the frequency domain about zero
for l=1:10 %Plot some realisations over non-negative frequencies
    plot(f(floor(n/2)+1:n),PSD(l,floor(n/2)+1:n),'LineWidth',0.1); 
end
plot(f(floor(n/2)+1:n),lowerCIPSD(floor(n/2)+1:n),'--r','LineWidth',2);
plot(f(floor(n/2)+1:n),upperCIPSD(floor(n/2)+1:n),'--r','LineWidth',2);
plot(f(floor(n/2)+1:n),PSDavg(floor(n/2)+1:n),'g','LineWidth',2.5); %Plots over non-negative frequencies, lines are too thick so can remove samples/remove thickness/reduce domain size
hold off
xlabel('Frequency');
ylabel('Power');
title('Power Spectral Density (Individual Realisations)');
xlim([0 1]);
subplot(2,1,2);
PSD2=abs(fftshift(fft(ACFavg)));
PSDlowCI=abs(fftshift(fft(lowerCI)));
PSDupCI=abs(fftshift(fft(upperCI)));
hold on
plot(f(floor(n/2)+1:n),PSDlowCI(floor(n/2)+1:n),'--r');
plot(f(floor(n/2)+1:n),PSDupCI(floor(n/2)+1:n),'--r');
plot(f(floor(n/2)+1:n),PSD2(floor(n/2)+1:n),'-b');
hold off
xlabel('Frequency');
ylabel('Power');
title('Power Spectral Density (Fourier transform of ACFavg)');
xlim([0 1]);

%% Sectioned out because takes a long time
RunGT1=2; %UPDATE this after you run this section to any other number (for speed)
if RunGT1==1
   q=zeros(1,length(lag)-2);
   [TrueCov,X1Cov,X2Cov,Q]=YCovariance(lag,40.099,40.099,10.025,10.025,25,q,0.01);
   disp('Done');
else 
    q=Q;
    Covariances=zeros(8,length(lag));
    figure(3);
    plot(dist,ACFavg,'k-');
    hold on 
    plot(dist,lowerCI,'--r');
    plot(dist,upperCI,'--r');
    for theta=2:2:16
        %[TrueCov,X1Cov,X2Cov,~]=YCovariance(lag,40.099,2.5,16.3,theta,25,q,0.01);
        [ICov,X1Cov,X2Cov,q]=TestIndicatorCov(lag,40.099,theta,16.3,10,-0.75,-0.5,40.099,3.5,25,1);
        Covariances(theta/2,:)=X1Cov;
        %TrueCorr=TrueCov/TrueCov(1); %Divide TrueCov by the variance to get the correlation
        %TrueCorr(length(TrueCorr))=0; %The last element is hugely negative for some reason, tried for many lags and always does the same thing 
        %plot(dist,TrueCorr);
        ICorr=ICov/ICov(1);
        ICorr(length(ICorr))=0;
        plot(dist,ICorr);
    end
    hold off
    ylim([-0.3 1.05]);
    xlim([0 70]);
    xlabel('Distance (mm)');
    ylabel('Autocorrelation');
    title('Comparison of data with model');
    legend('Average sample acf (data)','Lower limit of CI','Upper limit of CI','Theta=2','Theta=4','Theta=6','Theta=8','Theta=10','Theta=12','Theta=14','Theta=16','Location','eastoutside');
    figure(4);
    plot(dist,X2Cov,'-c');
    hold on
    for i=1:8
        plot(dist,Covariances(i,:));
    end
    hold off
    xlabel('Distance (mm)');
    ylabel('Autocorrelation');
    title('ACF of X1 and X2 processes');
    legend('X2','Theta=2','Theta=4','Theta=6','Theta=8','Theta=10','Theta=12','Theta=14','Theta=16','Location','eastoutside');
end

%ICov=TestIndicatorCov(lag,40.099,10.025,25);
%ICorr=ICov/ICov(1);
%ICorr(length(ICorr))=0;
%plot(dist, ICorr,'c-');

