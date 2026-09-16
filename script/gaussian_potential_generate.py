import numpy as np
import sys
data = np.loadtxt(sys.argv[1], dtype={'names': ('index', 'D_i', 'kuhn'), 'formats': ('U10', 'f8', 'f8')})
index = data['index']
D_i = data['D_i']
kuhn = data['kuhn']
Nisc = int(len(D_i))
npairs = sum(range(1, Nisc+1))
D_i= np.array(D_i)
D_pair_list = np.zeros((len(D_i)+1, len(D_i)+1))
for i in range(len(D_i)+1):
    for j in range(i, len(D_i)):
        D_pair_list[i][j] = (D_i[i]+D_i[j])/2


with open("input", 'w') as file:
    file.write("# Coarse-grained bottlebrush\nunits  lj\natom_style  angle\nboundary    p p p\n\n")
    file.write("# Minimization run\n\nread_data final0.cfg\nmass  * 1.0\n\npair_style gauss 40\n")
    for n in range(Nisc+1):
        for m in range(n, Nisc):
            coeff = D_pair_list[n][m]/2.278
            coeff = coeff**(-2)
            file.write("pair_coeff %d %d -3.12 %lf\n"%(n+1, m+1, coeff))
            # print("pair_coeff %d %d -3.12 %lf\n"%(n+1, m+1, coeff))
    file.write("\n")
    file.write("bond_style  harmonic\n")
    for n in range(Nisc-1):
        file.write("bond_coeff %d 200 %.2lf\n"%(n+1, (D_i[n]+D_i[n+1])/4))
    file.write("\n")
    file.write("angle_style harmonic\n")
    for n in range(Nisc-2):
        file.write("angle_coeff %d %lf 180\n"%(n+1,((kuhn[n]+kuhn[n+1])/((D_i[n]+D_i[n+1])/2)+(kuhn[n+2]+kuhn[n+1])/((D_i[n+2]+D_i[n+1])/2))/2))
    file.write("\n")
    file.write("neigh_modify every 1 delay 0 check yes\ntimestep    0.005\nthermo_style    custom etotal ebond temp ke eangle evdwl pe pxx press pyy vol pzz density\n")
    file.write("thermo_modify   line multi\n\n#min_style   cg\n#minimize 1.0e-8 1.0e-10 100000 1000000\n\n")
    file.write("# Main run\nfix   1 all nve\nfix   2 all langevin 1.0 1.0 1.0 157812\n\nthermo 10000\ndump  1 all custom 10000 dump.lammpstrj id mol type x y z\n")
    file.write("run   10000000\n\nwrite_data      final1.cfg\nwrite_restart   bb.restart\nclear")
    
    # write qsub.file
    with open("job.qsub", "w") as file:
        file.write("#!/usr/bin/env bash\n#$ -S /bin/bash\n#$ -j y\n")
        # file.write("#$ -q parallel2.q\n")
        file.write("#$ -q infiniband.q\n")
        file.write("#$ -cwd\n#$ -p -100\n#$ -pe mpi 4\n")
        file.write("#$ -N conn1\n")
        file.write("LMP_DIR=$SGE_O_HOME/app/lammps-stable_29Sep2021_update2/src\nmodule load gcc/9.3.0\nmpirun -np $NSLOTS $LMP_DIR/lmp_mpi -in input")
    