clc; clear all; close all;

[G, rock, p, testcase] = makeProblem(); 


fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);


inp = fullfile(fn, 'input');
info = dlmread(fullfile(inp, 'info.txt'));
    
Nf = info(1, 1);
Nc = info(2, 1);
offsets = info(3, :) + 1;
fine = dlmread(fullfile(inp, 'support.txt')) + 1;
coarse = rldecode((1:Nc), diff(offsets), 2);


