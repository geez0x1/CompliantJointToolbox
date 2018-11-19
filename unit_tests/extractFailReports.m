function [reports] = extractFailReports(testResults)
%EXTRACTFAILREPORTS Summary of this function goes here
%   Detailed explanation goes here

reports = {};
nRes = numel(testResults);

for iRes = 1:nRes
    if testResults(iRes).Failed ~= 0
%        reports{end+1} = testResults(iRes).Details.DiagnosticRecord.Report(:);
       reports{end+1} =  char(testResults(iRes).Details.DiagnosticRecord.EventLocation, testResults(iRes).Details.DiagnosticRecord.TestDiagnosticResults.DiagnosticText);
    end
    
end

celldisp(reports);

end

