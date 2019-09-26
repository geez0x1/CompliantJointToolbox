setCJTPaths;
outDir = ['.',filesep];
outFormat = 'pdf';
pubOpts = struct('format',outFormat,'outputDir',outDir);
publish('Take_A_Tour.m',pubOpts);