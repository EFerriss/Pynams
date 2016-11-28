% Ferriss 1/2015
% Saves peakfit information from FTIR_peakfit_loop
% to fname-peakfit.CSV

function savefit(output, spectraLocation, fname, file_ending, default_numPeaks)
% band positions, heights, widths, and areas
outputfile = [spectraLocation, fname, file_ending];
    fileID = fopen(outputfile, 'w');
formatSpec = '%i,%5.2f,%5.2f,%5.2f\n';
for row = 1:default_numPeaks
    fprintf(fileID, formatSpec, output(row,:));
end
fclose(fileID);
disp('fit saved')