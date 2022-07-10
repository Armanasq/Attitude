%% Main Script
% This script will call or attitude determination/estimation algorithms
% Algorithems:
% Attitude from Acc/Mag             AccMag
% Attitude from Acc/MagCal          AccMagCal

% Acc/Mag CF gyro (integral)        CFgyroIntg
% Acc/Mag CF gyro (Acc/Mag)         CFgyroAM
% Acc/Mag CF gyro (CF(t-1))         CFgyroCF
% TRIAD                             TRIAD
% Davenportds q-Method              davenport
% QUaternion ESTimator              QUEST
% Factored Quaternion Algorithm     FQA



clc
close all
clear all

%% Load Data
%load('slow_v4.mat');
load('BROAD_SampleRate.mat');
load('CleanDateNoNaN.mat');

fs = SampleRate;
input =input13;
output = output13;