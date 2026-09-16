import numpy as np
import sys 
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os

# Load files variables
r, pmf, mf = np.loadtxt(sys.argv[1], unpack=True) # reference potential
nsc1 = 10
d1 = nsc1**0.5*3.578
nsc_2 = 30
d2 = nsc_2**0.5*3.578
do = float(sys.argv[2]) # reference diameter
ro = r/do
lmdalist = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# AA interaction
nbin = 600
rcut = 0.03269*nbin
f = interp1d(ro*d1, pmf, kind='linear', fill_value='extrapolate')
g = interp1d(ro*d1, mf, kind='linear', fill_value='extrapolate')
new_r = np.arange(0.03269, rcut, 0.03269)
new_pmf = f(new_r)
new_mf = g(new_r)      
for n in range(len(new_pmf)):
    scaling = 1.4698682170542636
    # if new_pmf[n] > 0 and new_r[n] > (scaling*D_i[i]):
    if new_r[n] > (scaling*d1):
        new_pmf[n] = 0
        new_mf[n] = 0
r_AA = new_r
pmf_AA = new_pmf
mf_AA = new_mf

# BB interaction
f = interp1d(ro*d2, pmf, kind='linear', fill_value='extrapolate')
g = interp1d(ro*d2, mf, kind='linear', fill_value='extrapolate')
new_r = np.arange(0.03269, rcut, 0.03269)
new_pmf = f(new_r)
new_mf = g(new_r)
for n in range(len(new_pmf)):
    scaling = 1.4698682170542636
    # if new_pmf[n] > 0 and new_r[n] > (scaling*D_i[i]):
    if new_r[n] > (scaling*d2):
        new_pmf[n] = 0
        new_mf[n] = 0
r_BB = new_r
pmf_BB = new_pmf
mf_BB = new_mf
plt.plot(r_AA, pmf_AA, color="black", label='AA')
plt.plot(r_AA, pmf_BB, color="black", label='BB')

## MM/AM interaction - sigma scaling by lambda
pmfMM_lmda_list = []
mfMM_lmda_list = []
pmfAM_lmda_list = []
mfAM_lmda_list = []
for i in range(len(lmdalist)):
    # MM
    d_2MM = (1-lmdalist[i])*d1*d1+(lmdalist[i])*d2*d2
    d_MM = np.sqrt(d_2MM)
    new_r = np.arange(0.03269, rcut, 0.03269)
    f = interp1d(ro*d_MM, pmf, kind='linear', fill_value='extrapolate')
    g = interp1d(ro*d_MM, mf, kind='linear', fill_value='extrapolate')
    new_pmf = f(new_r)
    new_mf = g(new_r)
    for n in range(len(new_pmf)):
        scaling = 1.4698682170542636
        # if new_pmf[n] > 0 and new_r[n] > (scaling*D_i[i]):
        if new_r[n] > (scaling*d_MM):
            new_pmf[n] = 0
            new_mf[n] = 0
    r_MM = new_r
    pmf_MM = new_pmf
    mf_MM = new_mf
    pmfMM_lmda_list.append(pmf_MM)
    mfMM_lmda_list.append(mf_MM)
    # AM
    # d_2AM = 0.5*(d_2MM+d1*d1)
    # d_AM = np.sqrt(d_2AM)
    # new_r = np.arange(0.03269, rcut, 0.03269)
    # f = interp1d(ro*d_AM, pmf, kind='linear', fill_value='extrapolate')
    # g = interp1d(ro*d_AM, mf, kind='linear', fill_value='extrapolate')
    # new_pmf = f(new_r)
    # new_mf = g(new_r)
    # for n in range(len(new_pmf)):
    #     scaling = 1.4698682170542636
    #     # if new_pmf[n] > 0 and new_r[n] > (scaling*D_i[i]):
    #     if new_r[n] > (scaling*d_MM):
    #         new_pmf[n] = 0
    #         new_mf[n] = 0
    # r_AM = new_r
    pmf_AM = 0.5*(new_pmf+pmf_AA)
    mf_AM = 0.5*(new_mf+mf_AA)
    pmfAM_lmda_list.append(pmf_AM)
    mfAM_lmda_list.append(mf_AM)
    plt.plot(r_MM, pmf_MM, label='MM %lf'%(lmdalist[i]))
    plt.plot(r_MM, pmf_AM, label='AM %lf'%(lmdalist[i]))
    # plt.plot(r_MM, mf_MM, label='mfMM %lf'%(lmdalist[i]))
    # plt.plot(r_MM, mf_AM, label='mfAM %lf'%(lmdalist[i]))
    
plt.xlabel('r')
plt.ylabel('u(r)')
plt.legend()
plt.savefig("pot.pdf")
plt.show()

path = "/Users/haisukang/remote3/Workspace/SCMF_BB/shape/thermI_test_v1/f2nsc10-f2nsc30/v2_system_fix"
for l in range(len(lmdalist)):
    dir_main = os.path.join(path, f"l{lmdalist[l]:.1f}")
    dir_bb = os.path.join(dir_main, "bb")
    os.makedirs(dir_main, exist_ok=True)
    os.makedirs(dir_bb, exist_ok=True)

    # Write the pot0.table
    new_r = r_AA
    with open(os.path.join(dir_main, "pot0.table"), 'w') as file:
        file.write("# Bottlebrush melt ISC Vo for IBI pair AA, BB, AB, MM, AM\n\n")
        file.write("AA\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmf_AA[i]:.04f} {mf_AA[i]:0.04f}\n")
        file.write("\n")
        file.write("BB\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmf_BB[i]:.04f} {mf_BB[i]:0.04f}\n")
        file.write("\n")
        file.write("MM\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmfMM_lmda_list[l][i]:.04f} {mfMM_lmda_list[l][i]:0.04f}\n")
        file.write("\n")
        file.write("AM\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmfAM_lmda_list[l][i]:.04f} {mfAM_lmda_list[l][i]:0.04f}\n")
    with open(os.path.join(dir_bb, "pot0.table"), 'w') as file:
        file.write("# Bottlebrush melt ISC Vo for IBI pair AA, BB, AB, MM, AM\n\n")
        file.write("AA\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmf_AA[i]:.04f} {mf_AA[i]:0.04f}\n")
        file.write("\n")
        file.write("BB\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmf_BB[i]:.04f} {mf_BB[i]:0.04f}\n")
        file.write("\n")
        file.write("MM\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmfMM_lmda_list[l][i]:.04f} {mfMM_lmda_list[l][i]:0.04f}\n")
        file.write("\n")
        file.write("AM\nN %d R 0.001 %lf\n\n"%(len(new_r), rcut))
        for i in range(len(new_r)):
            file.write(f"{i+1} {new_r[i]:.06f} {pmfAM_lmda_list[l][i]:.04f} {mfAM_lmda_list[l][i]:0.04f}\n")

    # Write input file
    with open(os.path.join(dir_main,"input"), 'w') as file:
        file.write("# Coarse-grained bottlebrush\nunits  lj\natom_style  angle\nboundary    p p p\n\n")
        file.write("# Minimization run\n\nread_data final0.cfg\nmass  * 1.0\n\n")
        file.write("group M type 2\n")
        file.write("group A type 1\n")
        file.write("group Both intersect M A\n")
        file.write("group Monly subtract M Both\n")
        file.write("group Aonly subtract A Both\n\n")
        file.write("pair_style table linear %d\n"%(len(new_r)))  
        file.write("pair_coeff 1 1 pot0.table AA\n")           
        file.write("pair_coeff 2 2 pot0.table MM\n")
        file.write("pair_coeff 1 2 pot0.table AM\n")
        file.write("\n")
        file.write("bond_style  harmonic\n")
        file.write("bond_coeff 1 200.0 %lf\n"%(d1/2)) # for homopolymer systems
        blength_MM = (1-lmdalist[l])*d1/2 + lmdalist[l]*d2/2
        file.write("bond_coeff 2 200.0 %lf\n"%(blength_MM)) # for blend systems
        file.write("\n")
        file.write("angle_style harmonic\n")
        kuhn_1 = 0.067*2*d1
        kuhn_2 = 0.033*2*d2
        file.write("angle_coeff 1 %lf 180\n"%(kuhn_1)) # for homopolymer systems
        kuhn_MM = (1-lmdalist[l])*kuhn_1 + lmdalist[l]*kuhn_2
        file.write("angle_coeff 2 %lf 180\n"%(kuhn_MM)) # for blend systems
        file.write("\n")
        file.write("neighbor        2.0 bin\nneigh_modify    every 1 delay 0 check yes one 20000 page 2000000\ncomm_modify     cutoff 0.0   # >= your largest pair cutoff\n\n")
        file.write("neigh_modify every 1 delay 0 check yes\ntimestep    0.005\nthermo_style    custom etotal ebond temp ke eangle evdwl pe pxx press pyy vol pzz density\n")
        file.write("thermo_modify   line multi\n\nmin_style   cg\nminimize 1.0e-8 1.0e-10 100000 1000000\n\n")
        file.write("# Main run\nfix   1 all nve\nfix   2 all langevin 1.0 1.0 1.0 157812\n\n")

        file.write("compute pea all pe/atom bond angle\n")
        file.write("compute E_intra_M Monly reduce sum c_pea\n")  
        file.write("compute E_selfMM Monly group/group Monly pair yes\n")
        file.write("compute E_crossAM Monly group/group Aonly pair yes\n")
        file.write("variable Eintra  equal c_E_intra_M\n")
        file.write("variable EselfMM equal c_E_selfMM\n")
        file.write("variable Ecross  equal c_E_crossAM\n")
        file.write("thermo_style custom step temp v_Eintra v_EselfMM v_Ecross\n")
        file.write("thermo 10000\ndump  1 all custom 10000 dump.lammpstrj id mol type x y z\n")
        file.write("run   7000000\n\nwrite_data      final1.cfg\nwrite_restart   bb.restart\nclear")  
    with open(os.path.join(dir_bb,"input"), 'w') as file:
        file.write("# Coarse-grained bottlebrush\nunits  lj\natom_style  angle\nboundary    p p p\n\n")
        file.write("# Minimization run\n\nread_data final0.cfg\nmass  * 1.0\n\n")
        file.write("group M type 2\n")
        file.write("group A type 1\n")
        file.write("group Both intersect M A\n")
        file.write("group Monly subtract M Both\n")
        file.write("group Aonly subtract A Both\n\n")
        file.write("pair_style table linear %d\n"%(len(new_r)))  
        file.write("pair_coeff 1 1 pot0.table MM\n")
        file.write("pair_coeff 2 2 pot0.table MM\n")
        file.write("pair_coeff 1 2 pot0.table MM\n")
        file.write("\n")
        file.write("bond_style  harmonic\n")
        blength_MM = (1-lmdalist[l])*d1/2 + lmdalist[l]*d2/2
        file.write("bond_coeff 1 200.0 %lf\n"%(blength_MM)) 
        file.write("bond_coeff 2 200.0 %lf\n"%(blength_MM)) 
        file.write("\n")
        file.write("angle_style harmonic\n")
        kuhn_1 = 0.067*2*d1
        kuhn_2 = 0.033*2*d2
        kuhn_MM = (1-lmdalist[l])*kuhn_1 + lmdalist[l]*kuhn_2
        file.write("angle_coeff 1 %lf 180\n"%(kuhn_MM)) 
        file.write("angle_coeff 2 %lf 180\n"%(kuhn_MM)) 
        file.write("\n")
        file.write("neighbor        2.0 bin\nneigh_modify    every 1 delay 0 check yes one 20000 page 2000000\ncomm_modify     cutoff 0.0   # >= your largest pair cutoff\n\n")
        file.write("neigh_modify every 1 delay 0 check yes\ntimestep    0.005\nthermo_style    custom etotal ebond temp ke eangle evdwl pe pxx press pyy vol pzz density\n")
        file.write("thermo_modify   line multi\n\nmin_style   cg\nminimize 1.0e-8 1.0e-10 100000 1000000\n\n")
        file.write("# Main run\nfix   1 all nve\nfix   2 all langevin 1.0 1.0 1.0 157812\n\n")

        file.write("compute pea all pe/atom bond angle\n")
        file.write("compute E_intra_M Monly reduce sum c_pea\n")  
        file.write("compute E_selfMM Monly group/group Monly pair yes\n")
        file.write("compute E_crossAM Monly group/group Aonly pair yes\n")
        file.write("variable Eintra  equal c_E_intra_M\n")
        file.write("variable EselfMM equal c_E_selfMM\n")
        file.write("variable Ecross  equal c_E_crossAM\n")
        file.write("thermo_style custom step temp v_Eintra v_EselfMM v_Ecross\n")
        file.write("thermo 10000\ndump  1 all custom 10000 dump.lammpstrj id mol type x y z\n")
        file.write("run   7000000\n\nwrite_data      final1.cfg\nwrite_restart   bb.restart\nclear")
        
    # write qsub.file
    with open(os.path.join(dir_main,"job.qsub"), "w") as file:
        # Parallel
        file.write("#!/usr/bin/env bash\n#$ -S /bin/bash\n#$ -j y\n")
        file.write("#$ -q parallel.q\n")
        file.write("#$ -cwd\n#$ -p -100\n#$ -pe mpi 4\n")
        file.write("#$ -N TI_ISC%lf\n"%(lmdalist[l]))
        file.write("LMP_DIR=$SGE_O_HOME/app/lammps-stable_29Sep2021_update2/src\nmodule load gcc/9.3.0\nmpirun -np $NSLOTS $LMP_DIR/lmp_mpi -in input")
        # Serial
        # file.write("#$ -q xeon2.q,xeon1.q\n")
        # file.write("#$ -cwd\n#$ -p -100\n")
        # file.write("LMP_DIR=$SGE_O_HOME/app/lammps-stable_29Sep2021_update2/src\nmodule load gcc/9.3.0\n$LMP_DIR/lmp_serial -in input")
    with open(os.path.join(dir_bb,"job.qsub"), "w") as file:
        # Parallel
        file.write("#!/usr/bin/env bash\n#$ -S /bin/bash\n#$ -j y\n")
        file.write("#$ -q parallel.q\n")
        file.write("#$ -cwd\n#$ -p -100\n#$ -pe mpi 4\n")
        file.write("#$ -N TI_HISC%lf\n"%(lmdalist[l]))
        file.write("LMP_DIR=$SGE_O_HOME/app/lammps-stable_29Sep2021_update2/src\nmodule load gcc/9.3.0\nmpirun -np $NSLOTS $LMP_DIR/lmp_mpi -in input")
