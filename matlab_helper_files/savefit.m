% Ferriss 1/2015
% Saves peakfit information from FTIR_peakfit_afterpython
% to fname-peakfit.CSV

% band positions, heights, widths, and areas
outputfile = [spectraLocation, fname, '-peakfit.CSV'];
    fileID = fopen(outputfile, 'w');
formatSpec = '%i,%5.2f,%5.2f,%5.2f\n';
for row = 1:default_numPeaks
    fprintf(fileID, formatSpec, output(row,:));
end
fclose(fileID);
