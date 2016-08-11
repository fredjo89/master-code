
% Delete previous mex file
clc; clear all; close all;
aString = [pwd,'/mex_cppMultiscaleBasis.mexa64'];
delete(aString);
disp('Mex Deleted'); 
aString = [pwd,'/myWork/my_mex_cppMultiscaleBasis.mexa64'];
delete(aString);
disp('Mex Deleted'); 
