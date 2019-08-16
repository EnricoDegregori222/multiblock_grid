clear all
clc
close all

mbk=7;
% Grid
lxi0 = 230;
lxi01 = 140;
lxi1 = 800;
lxi2 = 400;
let0 = 400;
let01 = 300;
let1 = 50;
lze0 = 1;
smg = 4e-5; % Re = 500000
% smg = 1e-4; % Re = 200000
smgvr = 1.0;
fracvr = 3/20;
ptvr = 3/5;
spx = 0.0;
spy = 0.0;

lxie0=lxi0;
lxis1=lxie0+1;
lxie1=lxis1+lxi1;
lxis2=lxie1+1;
lxit = lxi0+lxi1+lxi2+2;

lete0=let0;
lets1=lete0+1;
lete1=lets1+let1;
lets2=lete1+1;
lett=2*let0+1;
lettr=2*let0+let1+2;

lete01=let01;
lets11=lete01+1;
lete11=lets1+let1;
lets21=lete11+1;
lett1=2*let01+1;
lettr1=2*let01+let1+2;

% Domain
domh = 12.0;
domh1 = 1.0; % 1.5 for straight LE / 1.0 for wavy LE
skew = 0.0;
doml0 = 7.0;
doml01 = 0.5;
doml1 = 12.0;
span = 0.05;

% WLE
wlew = 0.05;
wlea = 0.0;

% Sponge
szth0 = 2.0;
szth1 = 2.0;
nsponge = 16;

% Constants
free = 1e6;

% Airfoil section mesh
filename_input = 'v2c_ref.dat';
nref = 3;
n_initial_points = 8;

% Grid save
gridsave = 0;
gridtec = 1;