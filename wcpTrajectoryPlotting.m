dataDir = 'F:\WCP_Traj';
fn1 = '20250909_GLSAR_Lift1-FIELD.txt';
fn2 = '20250909_GLSAR_Lift2-FIELD.csv';
fn3 = '20250909_GLSAR_Lift3-FIELD.csv';
fn4 = '20250909_GSSAR_Lift1-FIELD.csv';
fn5 = '20250909_GSSAR_Lift2-FIELD.csv';

dataDir2 = 'F:\WCP_Traj\new_GSSAR';

traj1 = readtable([dataDir,'\',fn1]);
traj2 = readtable([dataDir,'\',fn2]);
traj3 = readtable([dataDir,'\',fn3]);
traj4 = readtable([dataDir,'\',fn4]);
traj5 = readtable([dataDir,'\',fn5]);

% new GSSAR
fn6 = '20250909_GSSAR_Lift1-FIELD.txt';
fn7 = '20250909_GSSAR_Lift2-FIELD.txt';
traj6 = readtable([dataDir,'\',fn6]);
traj7 = readtable([dataDir,'\',fn7]);

% rtk GGSAR
fn8 = '20250908_GSSAR_Lift1-FIELD.txt';
fn9 = '20250908_GSSAR_Lift2-FIELD.txt';
traj8 = readtable([dataDir,'\',fn8]);
traj9 = readtable([dataDir,'\',fn9]);

% rtk GLSAR
fn10 = '20250908_GLSAR_Lift1-FIELD.txt';
fn11 = '20250908_GLSAR_Lift2-FIELD.txt';
traj10 = readtable([dataDir,'\',fn10]);
traj11 = readtable([dataDir,'\',fn11]);

% 091025
dataDir3 = 'F:\WCP_Traj\20250910';
fn12 = '20250910_north_GLSAR_Lift1-FIELD.txt';
fn13 = '20250910_north_GLSAR_Lift2-FIELD.txt';
fn14 = '20250910_south_GLSAR_Lift1-FIELD.txt';
fn15 = '20250910_south_GLSAR_Lift2-FIELD.txt';
fn16 = '20250910_south_GSSAR_Lift1-FIELD.txt';
fn17 = '20250910_south_GSSAR_Lift2-FIELD.txt';
fn18 = '20250910_north_GSSAR_Lift1-FIELD.txt';
fn19 = '20250910_north_GSSAR_Lift2-FIELD.txt';

traj12 = readtable([dataDir3,'\',fn12]);
traj13 = readtable([dataDir3,'\',fn13]);
traj14 = readtable([dataDir3,'\',fn14]);
traj15 = readtable([dataDir3,'\',fn15]);
traj16 = readtable([dataDir3,'\',fn16]);
traj17 = readtable([dataDir3,'\',fn17]);
traj18 = readtable([dataDir3,'\',fn18]);
traj19 = readtable([dataDir3,'\',fn19]);
