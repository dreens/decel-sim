load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/ModeSeqSlow142_1_October-27-2015_17-40-41/results_ModeSeqSlow142_1.mat')
r1 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/ModeSeqSlow142_2_October-27-2015_17-40-42/results_ModeSeqSlow142_2.mat')
r2 = r;
r = [r1 r2];
clear r1 r2
results_initfinalphase(r)



load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/ModeSeqTrap142_1_October-27-2015_20-52-55/results_ModeSeqTrap142_1.mat')
r1 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/ModeSeqTrap142_2_October-27-2015_20-52-57/results_ModeSeqTrap142_2.mat')
r2 = r;
r = [r1 r2];
clear r1 r2
any(r(1).modeseq==3)
r(1).molnum(end)
any(r(2).modeseq==3)
r(2).molnum(end)
results_synch(r)


load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/ModeSeqT440_2_October-27-2015_21-37-05/results_ModeSeqT440_2.mat')
r1 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/ModeSeqT440_1_October-27-2015_21-37-03/results_ModeSeqT440_1.mat')
r2 = r;
r = [r1 r2];
clear r1 r2
any(r(1).modeseq==3)
r(1).molnum(end)
any(r(2).modeseq==3)
r(2).molnum(end)
results_synch(r)


load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_phase_Gaussian_1_October-30-2015_18-44-36/results_S=13comp_phase_Gaussian_1.mat')
r1 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_phase_Gaussian_2_October-30-2015_18-44-37/results_S=13comp_phase_Gaussian_2.mat')
r2 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_phase_Gaussian_3_October-30-2015_18-44-38/results_S=13comp_phase_Gaussian_3.mat')
r3 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_phase_Gaussian_4_October-30-2015_18-44-39/results_S=13comp_phase_Gaussian_4.mat')
r4 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_phase_Gaussian_5_October-30-2015_18-44-40/results_S=13comp_phase_Gaussian_5.mat')
r5 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_phase_Gaussian_6_October-30-2015_18-44-41/results_S=13comp_phase_Gaussian_6.mat')
r6 = r;
r = [r1 r2 r3 r4 r5 r6];
results_compareEricphase(r)



load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_bunch_1_October-31-2015_12-33-42/results_S=13comp_bunch_1.mat')
r1 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_bunch_2_October-31-2015_12-33-43/results_S=13comp_bunch_2.mat')
r2 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_bunch_3_October-31-2015_12-33-44/results_S=13comp_bunch_3.mat')
r3 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_bunch_4_October-31-2015_12-33-45/results_S=13comp_bunch_4.mat')
r4 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_bunch_5_October-31-2015_12-33-46/results_S=13comp_bunch_5.mat')
r5 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/S=13comp_bunch_6_October-31-2015_12-33-47/results_S=13comp_bunch_6.mat')
r6 = r;
r = [r1 r2 r3 r4 r5 r6];
results_compareEricphase(r)


%%
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/try_widening_1_November-02-2015_10-19-11/results_try_widening_1.mat')
r1 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/try_widening_2_November-02-2015_10-19-12/results_try_widening_2.mat')
r2 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/try_widening_3_November-02-2015_10-19-13/results_try_widening_3.mat')
r3 = r;
load('/users/ye/dare4983/Documents/MATLAB/slowANDtrap/Cluster/try_widening_4_November-02-2015_10-19-14/results_try_widening_4.mat')
r = [r1 r2 r3 r];
results_compareEricphase(r([1 3 2 4]))