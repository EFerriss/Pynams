% Ferriss 1/2015
% take a baseline-subtracted spectrum 
% generated using python FTIRspectra.py self.save_baseline()
% and fits the peaks using peakfit from the Matlab file exchange
% Use savefit.m to save fit to filename-peakfit.csv
%
clear all; close all; 

%% User input section 

% where the files are stored
spectraLocation = 'C://Users/Ferriss/Documents/Code/olivine/FTIR/';
saveFigureLocation = spectraLocation;

% ending of file where FTIR baseline is stored
baseline_ending = '-baseline-low.CSV';

% ending of file where peakfit information will be saved
file_ending = '-peakfit-low.CSV';

fname_list = {'SC1-2-800C-7hr-a01'}

%default_peakpos = [3645 3617 3540 3460 3443 3355 3305 3250]; % diopsides
%default_peakpos = [3620 3550 3460 3355]; % PMR-53
default_peakpos = [3572 3525 3482.8 3600. 3542 3566.0]; % San Carlos olivine h
%  3325 3355 

default_peakpos = sort(default_peakpos);

for r=1:length(fname_list)
    
    disp(' ')
    fname = char(fname_list(r))
    use_previous_fit = 0; % 1 not working right now

    set_peakpos = 0; % = 1 to set peak positions manually
        changePeak = 1:6;
        manual_peakpos = [3443 3355 3305 3250];
        %[3645 3617 3540 3460 3443 3355];
        
    set_heights = 0; % = 1 to set peak height manually
    set_widths = 0; % = 1 to set peak height manually

        changePeakH =  [3443 3355 3305 3250];
        manual_heights = [2.1 1. 0. 0.];
        
        changePeakW  = [3443 3355 3305 3250];
        manual_widths = [17. 45. 0. 0.];
        
        % 16 = fixed-position Gaussian
    % 17 = fixed-position Lorentzian
    peakshape = 16;
    
    set_title = 0;
        manual_title = '';
    
        
% -----------------------------------------------------------------
% End user input section
% -----------------------------------------------------------------
    pek_outputfile = [spectraLocation, fname, baseline_ending];
    if exist(pek_outputfile) == 0
       disp('You have to subtract off a baseline in python');
       disp('and save with self.save_baseline() to make');
       disp([fname, baseline_ending])
       
       return
    else
       previous = xlsread(pek_outputfile);
       wn = previous(2:end,1);
       % previous(2:end,2) are baseline absorbance values
       pek = previous(2:end,3); 
       wi = min(wn);
       wf = max(wn);
    end

    if use_previous_fit == 1
        peakfitfile = [spectraLocation, fname, file_ending];
       if exist(peakfitfile) == 0
           disp('Fit peaks, then use savefit.m to save fname-peakfit.CSV');
           return
       else
          disp('Using previous fit from fname-peakfit.CSV');       
          previous_fit = xlsread(peakfitfile);
          peakpos = previous_fit(:,1);
          numPeaks = length(peakpos);
          heights = previous_fit(:,2);
          widths = previous_fit(:,3);
          fitpeakarea = previous_fit(:,4);
          peakshape = previous_fit(1,5);
          % Bell calibration, ppm H2O
          disp('water estimate based on summed band areas');
          water = sum(fitpeakarea) * 3 * 0.141;
       end
       disp(water)
    else

        if set_peakpos == 1
            disp('Manually setting peak positions')
            peakpos = manual_peakpos;
        else
            disp('Using default peak positions')
            peakpos = default_peakpos;
        end

        numPeaks = length(peakpos);

        [DR,Derror]=...
                peakfit([wn pek],0,0,numPeaks,peakshape,0,0,0,0,peakpos);

        % fit results
        peakpos = DR(:,2);
        heights = DR(:,3);
        widths = DR(:,4);
        fitpeakarea = DR(:,5);

        if set_heights == 1
            disp('Manually setting heights')
                for k = 1:length(changePeakH)
                    index = find(peakpos == changePeakH(k));
                    if isempty(index)
                        disp('Check peak position values in changePeakH');
                    else
                        heights(index) = manual_heights(k);
                    end
                end
        end

        if set_widths == 1
            disp('Manually setting widths')
                for k = 1:length(changePeakW)
                    index = find(peakpos == changePeakW(k));
                    if isempty(index)
                        disp('Check peak position values in changePeakW');
                    else
                        widths(index) = manual_widths(k);
                    end
                end
        end

        % Set up output to include ALL default peak positions for consistency
        % across various profiles. Override values where peakpos = default_peakpos.
        default_numPeaks = length(default_peakpos);
        output = zeros(default_numPeaks, 4);
        output(:, 1) = default_peakpos;

        for k = 1:numPeaks
            % Clean up scraggly numbers
            thresh = 0.005;
            threshhigh = 400.;    
            if (heights(k) <= thresh || widths(k) <= thresh || ...
                fitpeakarea(k) <= thresh || widths(k) >= threshhigh || ...
                heights(k) >= threshhigh)
                heights(k) = 0.0;
                widths(k) = 0.0;
                fitpeakarea(k) = 0.0;
            end
                      
            index = find(default_peakpos == peakpos(k));
            output(index, 2) = heights(k);
            output(index, 3) = widths(k);
            output(index, 4) = fitpeakarea(k);
        end
    end % end if statement for using previous fit

%% Plotting details you might want to play with

% y-axis limits in the plot
bottom = 0;
top = (max(pek))+ 0.35*(max(pek));

% where to put peak wavenumber labels along y-axis
shift = top - top*0.2;

% jpeg resolution
resolution = '-r100';
numbersize = 18;
linesize = 2;
labelsize = 20;
legendtext = 16;

measuredline.color = [0.5 0.5 0.5];
measuredline.linewidth = 5;
measuredline.marker = 'none';
measuredline.linestyle = '-';

calculatedline.color = 'k';
calculatedline.linewidth = 2;
calculatedline.marker = 'none';
calculatedline.linestyle = '--';

bandline.color = [0 0.6 0];
bandline.marker = 'none';
bandline.linestyle = '-';
bandline.linewidth = 2;

accent.linestyle = '-';
accent.color = 'k';
accent.linewidth = 1;
      
jcolor(1,:) = [0 0 1];
jcolor(2,:) = [1 0 0];
jcolor(3,:) = [0 0.5 1];
jcolor(4,:) = [0 0.6 0];
jcolor(5,:) = [1 0.5 0];
jcolor(6,:) = [1 0 1];
jcolor(7,:) = [1 1 0];
jcolor(8,:) = [0 1 0.5];
jcolor(9,:) = [0.5 0 0];
jcolor(10,:) = [0.5 0 1];
jcolor(11,:) = [1 0 0];

labeltextformat.fontsize = legendtext;
labeltextformat.rotation = 90;
labeltextformat.backgroundcolor = 'w';
labeltextformat.horizontalalignment = 'left';
labeltextformat.color = 'k';

w = 900; h = 520; % figure pixel dimensions

% figure
close all;

% gcfpos = [25 100 w h];
%gcfpos = [25 500 900 h];
gcfpos = [25 35 900 h];

set(gca,'FontSize',numbersize,'Nextplot','add',...
        'Box','on','XGrid','off','YGrid','on',...
        'LineWidth',linesize,'Xdir','reverse');
set(gcf,'Color','w','position',gcfpos);

% dummies for legend
plot(7000,0,measuredline); 
plot(7000,0,calculatedline); 
plot(7000,0,bandline); 

% line labels
for j = 1:numPeaks
    if heights(j) ~= 0 && widths(j) ~= 0
    line([peakpos(j) peakpos(j)],[bottom top],accent,...
         'Color',jcolor(j,:));
    text(peakpos(j),bottom+shift,num2str(peakpos(j)),...
        labeltextformat);
    end
end

LorT = zeros(length(wn),1); % prep total calculated curve

clear yLor;

fitpeakarea_manual = zeros(numPeaks, 1);
for k = 1:numPeaks
    % determine individual Gaussians
    yLor(:,k)=heights(k).*exp(-((wn-peakpos(k))./...
              (0.6005615.*widths(k))).^2);
        
    % plot individual bands
    plot(wn,yLor(:,k),bandline,'Color',jcolor(k,:)); 

    fitpeakarea_manual(k) = sum(yLor(:,k)) * (wn(2)-wn(1));
    
    % Add onto total curve = sum of individual Lorenztians
    LorT=LorT+yLor(:,k);
    
end


% Changing the area in the output when manualy fiddling
% Note the length of the output, which is always based on the 
% default peak positions, may well differ from that of the
% actual peakpos peak positions and so requires a separate 
% indexing step to get the correct value from fitpeakarea_manual
if set_heights == 1
    for ktan = 1:length(changePeakH)
        index = find(output(:,1) == changePeakH(ktan));
        areaindex = find(peakpos == changePeakW(ktan));
        output(index, 4) = fitpeakarea_manual(areaindex);
    end
end

if set_widths == 1
    for ktan = 1:length(changePeakW)
        index = find(output(:,1) == changePeakW(ktan));
        areaindex = find(peakpos == changePeakW(ktan));
        output(index, 4) = fitpeakarea_manual(areaindex);
    end
end

disp('fit results:')
disp(output)

% bulk measured curve
plot(wn,pek,measuredline);
% summation curve
plot(wn,LorT,calculatedline); 

axis([wi wf bottom top]);

xlabel('wavenumber (cm^{-1})','Fontsize',labelsize);
ylabel('absorbance (cm^{-1})','Fontsize',labelsize);

if set_title == 1
    title(manual_title,'FontSize',labelsize);
else
    title(fname,'Fontsize',labelsize);
end
    
str1=sprintf('individual\nbands');
str2=sprintf('calculated\nspectrum');
hleg=legend('observed',str2,str1);
set(hleg,'LineWidth',linesize,'FontSize',legendtext,...
        'Location','Northeast');

%Fig=[saveFigureLocation, fname, file_ending, '.jpg'];
%export_fig(Fig, '-jpeg',resolution,'-painters');

savefit(output, spectraLocation, fname, file_ending, default_numPeaks)
    
end % end loop through all spectra in fname_list