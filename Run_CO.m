outDir = ['..',filesep, 'results'];
outFormat = 'pdf';
pubOpts = struct('format',outFormat,'outputDir',outDir);
publish('Take_A_Tour_CO.m',pubOpts);