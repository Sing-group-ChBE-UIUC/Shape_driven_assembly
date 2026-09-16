import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import root_scalar
from pathlib import Path


# Define Nsc(x) as a function — replace with your desired functional form
# #Tapered
# def Nsc_1(x):
    # return -0.6 * x/0.6 + 70  # Example from image
    # return -0.3015*x/0.6+70.3015

# def Nsc_1(x, nchange):
#     # return -0.6 * x/0.6 + 70  # Example from image
#     # return -0.3015*x/0.6+70.3015
#     # list of shape changing
#     # slope_list = np.linspace(-60/99/1.2, 0, nchange)
#     slope_list = [-0.02, -0.06, -0.09]
#     contour=[]
#     for i in range(nchange):
#         contour.append(slope_list[i]*(x-60)+45)
#     return contour

# # diamond - half and symmetry
# def Nsc_1(x):
#     # return 1.2245*x/0.6+8.7755
#     return 0.606*x/0.6+10-0.606
# def Nsc_2(x):
#     # return -1.2245*x/0.6+131.225
#     return -0.6*x/0.6+130
# def Nsc_1(x, nchange):
#     # slope_list = np.linspace(0, 60/49, nchange)
#     slope_list = np.linspace(0, 60/99/1.2, nchange)
#     slope_list = slope_list[::-1]
#     contour=[]
#     for i in range(nchange):
#         contour.append(slope_list[i]*(x-60)+40)
#     return contour
# def Nsc_2(x, nchange):
#     slope_list = np.linspace(-60/99/1.2,0, nchange)
#     contour=[]
#     for i in range(nchange):
#         contour.append(slope_list[i]*(x-180)+40)
#     return contour
# # #bowtie
# def Nsc_1(x, nchange):
#     # slope_list = np.linspace(-60/49/1.2,0, nchange)
#     slope_list = np.linspace(-60/99/1.2,0, nchange)
#     contour=[]
#     for i in range(nchange):
#         contour.append(slope_list[i]*(x-60)+40)
#     return contour
# def Nsc_2(x, nchange):
#     # slope_list = np.linspace(0,60/49/1.2, nchange)
#     slope_list = np.linspace(0,60/99/1.2, nchange)
#     slope_list = slope_list[::-1]
#     contour=[]
#     for i in range(nchange):
#         contour.append(slope_list[i]*(x-180)+40)
#     return contour
# #football
def Nsc_1(x):
    # return -0.0245*(x/0.6)**2+2.4735*x/0.6+7.551
    y_ini=-0.0245*(x/1.2)**2+2.4735*(x/1.2)+7.551
# hourglass
# def Nsc_1 (x):
# #    return -0.1*(x/0.6)**2+5.1*x/0.6+5
#    return -0.1*(x/1.2)**2+5.1*x/1.2+5
# def Nsc_2 (x):
# #    return -0.1*(x/0.6+1)**2+15.1*(x/0.6+1)-500
#    return -0.1*(x/1.2)**2+15.1*x/1.2-500
# #cylinder
# def Nsc_1(x):
#     return 70

# Parameters
nsc_min = 10
nsc_max = 70
alpha = 3.578  # Scaling constant
x_max = 120 # End of backbone #football
tolerance = 1e-2  # To stop adding beads once we're close to the end
nchange = 20
# slope_list = np.linspace(-0.5025/3*2, 0, nchange)  
# slope_list1 = np.linspace(0, 60/99/1.2, nchange)
# # slope_list1 = slope_list1[::-1]
slope_list1= np.linspace(-60/99/1.2, 0, nchange)
# # slope_list2= np.linspace(-60/49/1.2, 0, nchange)
# slope_list1 = slope_list1[::-1]
# slope_list1 = [-0.02, -0.06, -0.09]
# Other continuous contour 
x_start = 0 #bowtie
# x_start = -4.1#diamond
D_list = []
x_list = []
# Define the equation to solve: f(D) = D^(3/2) - alpha * integral
# def make_equation1(x_start, slope):
#     # Nsc(x) = slope*(x-30) + 45
#     # Integral can be done analytically too, but keep quad for now.
#     def f(D):
#         val, _ = quad(lambda x: slope*(x - 30.0) + 40.0, x_start, x_start + D)
#         return D**3 - alpha**2 * val
#     return f
# def make_equation2(x_start, slope):
#     # Nsc(x) = slope*(x-30) + 45
#     # Integral can be done analytically too, but keep quad for now.
#     def f(D):
#         val, _ = quad(lambda x: slope*(x - 90.0) + 40.0, x_start, x_start + D)
#         return D**3 - alpha**2 * val
#     return f

def solve_D_list_for_slope(slope, x_start0, x_max, alpha):
    def make_equation1(x_start):
        def f(D):
            val, _ = quad(lambda x: slope*(x - 60.0) + 40.0, x_start, x_start + D)
            return D**3 - alpha**2 * val
        return f

    x_start = float(x_start0)
    D_list = []
    x_list = []

    while x_start <= x_max:
        f = make_equation1(x_start)
        try:
            sol = root_scalar(f, bracket=[0.1, 100.0], method="brentq")
        except ValueError:
            break
        if sol.converged:
            if x_start + sol.root > 130:
                break
            else:
                D_i = sol.root
                # print(D_i)
                x_list.append(x_start)
                D_list.append(D_i)
                x_start += D_i
        else:
            break

    # end position (where the last segment ends)
    if len(D_list) == 0:
        x_end = x_start0
    else:
        x_end = x_list[-1] + D_list[-1]
    return np.array(x_list), np.array(D_list), x_end


# def equation1(D):
#     # if D < 0:
#     #     return 1e6  # Prevent negative or zero roots
#     # else:
#         integral1= quad(lambda x: Nsc_1(x, nchange)[nth], x_start, x_start+D)
#         return D**(3) - alpha**2 * (integral1[0])
# def equation2(D):
#     # if D < 0:
#     #     return 1e6  # Prevent negative or zero roots
#     # else:
#         integral1= quad(lambda x: Nsc_2(x, nchange)[nth], x_start, x_start+D)
#         return D**(3) - alpha**2 * (integral1[0])
# D_list_nchange = []
# x_list_nchange = []
# ISC_list_nchange = []

# print(slope_list1)
for nth in range(nchange):
    plt.figure(figsize=(4,7))
    ax=plt.subplot(1,1,1)
    slope = slope_list1[nth]
    x_start0 = 0

    # 1st pass
    x_list, D_list, x_end = solve_D_list_for_slope(slope, x_start0, x_max, alpha)

    # # If overshoot, adjust start and recompute
    for _ in range(6):
        if x_end > x_max and len(D_list) > 0:
            dev = x_end - x_max
            x_start0_adj = x_start0 - 3   # shift left so it ends at x_max
            x_list, D_list, x_end2 = solve_D_list_for_slope(slope, x_start0_adj, x_max, alpha)
            print(f"nth {nth}: overshoot {dev:.4f}, recomputed with x_start0={x_start0_adj:.4f}, old_end ={x_end:.2f}, new end={x_end2:.4f}, Nisc={len(D_list)}")
            x_end = x_end2
            x_start0 = x_start0_adj
        if x_end <= 120:
            break
        
        if x_end < x_max:
            dev = x_max - x_end
            x_start0_adj = x_start0 + dev   # shift left so it ends at x_max
            x_list, D_list, x_end2 = solve_D_list_for_slope(slope, x_start0_adj, x_max, alpha)
            # print(f"nth {nth}: undershoot {dev:.4f}, recomputed with x_start0={x_start0_adj:.4f}, new end={x_end2:.4f}")
            x_end = x_end2


    # while x_start < x_max:
    #     try:
    #         sol = root_scalar(equation2, bracket=[0.1, 100], method='brentq')
    #         if sol.converged:
    #             D_i = sol.root
    #             x_list.append(x_start)
    #             D_list.append(D_i)
    #             x_start += D_i
    #             print(D_i)
    #         else:
    #             break
    #     except ValueError:
    #         # If no root found in bracket range
    #         print("root no range")
    #         break
    # D_list = D_list + D_list[::-1]
    D_list = np.concatenate([D_list, D_list[::-1]])
    x_list_1 = []
    x_start = x_start0
    for i in range(len(D_list)):
        x_list_1.append(x_start)
        x_start += D_list[i]
    D_list = np.array(D_list)
    x_list = np.array(x_list)
    print(len(D_list))
    # add middle bead
    Tot_ISC = len(D_list)+len(D_list)-1
    ISC_list = []
    # homopolymer
    for i in range(1,Tot_ISC+1):
        if i%2 != 0:
            D_ISC = D_list[int((i)/2)]
            ISC_list.append(D_ISC)
        else: 
            D_ISC = (D_list[int((i-1)/2)]+D_list[int((i+1)/2)])/2
            ISC_list.append(D_ISC)
    mother_path = Path(f"/Users/haisukang/remote3/Workspace/SCMF_BB/shape/v3_evenNisc_longbb/football_70_10/shape_change/s{slope_list1[nth]:.2f}")
    mother_path.mkdir(parents=True, exist_ok=True)
    file_path = mother_path/f"D_cfg_{slope_list1[nth]:.2f}.txt"
    with open(file_path, 'w') as file:
        for i in range(Tot_ISC):
            file.write("isc%d %.2lf %.4lf\n"%(i+1, ISC_list[i], 2.2534*ISC_list[i]-6.7449))
               

    cmap1 = plt.get_cmap('Oranges')
    # Plotting result similar to your figure
    x_vals1 = np.linspace(0, 120, 200)
    x_vals2 = np.linspace(120, 240, 200)
    # x_vals1 = np.linspace(0, 60, 200)
    # x_vals2 = np.linspace(60, 125, 200)
    y_contour1 = Nsc_1(x_vals1, nchange)
    y_contour2 = Nsc_2(x_vals2, nchange)
    y1 = y_contour1[nth]
    y2 = y_contour2[nth]

    ax.plot(x_vals1, y1, '--', color = 'black', label=r'$N_{sc}(x)$')
    ax.plot(x_vals2, y2, '--', color = 'black')
    ax.plot(x_vals1, 3.578*y1**0.5, '-', color = 'purple', label=r'$\alpha N_{sc}^{1/2}$')
    ax.plot(x_vals2, 3.578*y2**0.5, '-', color = 'purple')

    for i, (x, D) in enumerate(zip(x_list_1, D_list)):
        circle = plt.Circle((x+D/2, D/2), D/2, color='lightblue', ec='purple', fill=False)
        ax.add_patch(circle)

    # ax.set_xlim(0, x_max + 10)
    ax.set_ylim(-5, 80) 
    ax.set_xlabel(r"$N_{bb}$")
    ax.set_ylabel(r"$N_{sc}$ or D")
    ax.set_aspect('equal')
    ax.legend(ncol=3)

    plt.savefig(mother_path/f"ISC_s{slope_list1[nth]:.2f}.pdf", bbox_inches="tight")
    plt.show()
    plt.close()

