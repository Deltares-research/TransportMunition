# Transport model for UXO

# Packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

def transport(params, files):


    # Input
    pathfile = ".\\Excel\\"
    outputpath = ".\\Plaatjes\\"
    outputxcl = pathfile+files["output_file"]
    rho_w = params["rho_w"]
    g = params["g"]
    mu = params["mu"]
    factr =  1
    pp = 1
    dbedtok =1
    percent_bound = 0.1 # 20% bound
    threshold =1 / 50
    dragCoeff = 1.75

    # OOlijst
    dt = pd.read_excel(pathfile+files["OOlijst"])
    dt = dt.to_numpy()
    dt_name = dt[:, -3:-1]
    dt = dt[:, :-3].astype('float')

    dk = pd.read_excel(pathfile+files["ruwheid"], header=None)
    dk = dk.to_numpy()
    header = dk[0, :]
    indx = np.unique(header, return_index=True)
    indx = indx[1][:-1]
    objects = dk[1:, 0]
    dk = dk[1:,indx]
    no = np.size(dk, 1)
    header = header[indx].astype("str")
    header = np.char.replace(header, " [m]", "")
    minSnel = []
    minSnel.append(objects.astype("str"))
    column = []
    column.append("Objects")
    for ii in range(0, no):
        # Object shape for RWS
        Dmax = dt[:, 0]
        L = dt[:, 1]
        V_obj = dt[:, 4]
        rho_obj = dt[:, 3]
        k = dk[:, ii].astype("float")

        # Formulas
        "Menzel"
        a1 = 0.2315 # Menzel
        b1 = 2.1137 # Menzel
        D_avg = np.sqrt(4*V_obj/(np.pi*L))
        k[k==0] = 2.32*1e-6  # At least have non zero
        zbs = k/Dmax
        U_crit = ((rho_w*9.81*D_avg*V_obj*(rho_obj-rho_w))/(2*a1*mu**2*L))**(1/b1)*mu/(rho_w*(1-zbs)*D_avg)
        thetatUcrit = U_crit**2/(g*D_avg*(rho_obj/rho_w-1)) *(D_avg-k) / D_avg

        "Rennie"
        a2 = 1
        b2 = -0.75
        alpha2 = 25.3*(D_avg/(factr*k))**b2
        thetatUcrit2 = 3.3 * a2 * np.tan(alpha2/180*np.pi)/(1+0.8*np.tan(alpha2/180*np.pi))
        U_crit2 = np.sqrt(thetatUcrit2*(g*D_avg*(rho_obj/rho_w-1))*D_avg/(D_avg-factr*k))

        "Friedrichs"
        thetatUcrit3 = 1.2*(D_avg/(factr*k))**(-0.72)
        U_crit3 = np.sqrt(thetatUcrit3*(g*D_avg*(rho_obj/rho_w-1))*D_avg/(D_avg-factr*k)) # Sediment transport: 2016a Friedrichs

        "Rennie"
        thetatUcrit4 = 1.64*(D_avg/(factr*k))**(-0.71)
        U_crit4 = np.sqrt(thetatUcrit4*(g*D_avg*(rho_obj/rho_w-1)) *D_avg/(D_avg-factr*k))

        "Chu et al - seems more comparable to Menzel"
        b=np.ones_like(k) *factr* k*0.4445 # based on Chu et al. 2022, around factor 2 smaller than alpha2, is angle
        gamma =0.2
        Cd =2 # Conservative value.
        term = ((1-b/D_avg)/np.sqrt(b/D_avg*(1-b/D_avg)))*(1+factr*k/D_avg-2*b/D_avg)+2*gamma
        U_crit5 = np.sqrt(g * D_avg*(rho_obj/rho_w-1)/Cd*np.pi/term)
        thetatUcrit5 = U_crit5**2/(g*D_avg*(rho_obj/rho_w-1)) *(D_avg-factr*k) / D_avg

        "Martin"
        a1 = 0.25
        b1 = 2 # Found on old database
        D_avg = np.sqrt(4*V_obj/(np.pi*L))
        k[k==0] = 2.32*1e-6  # At least have non zero
        zbs = k/Dmax
        # dragCoeff = 8.1744*D_avg/Dmax -5.6744
        U_crit6 = ((rho_w*9.81*Dmax*V_obj*(rho_obj-rho_w))/(2*a1*mu**2*L))**(1/b1)*mu/(rho_w*(1-zbs)*Dmax) / dragCoeff
        thetatUcrit = U_crit**2/(g*Dmax*(rho_obj/rho_w-1)) *(Dmax-k) / Dmax

        # Gathering data
        U_crit_tot=np.array([U_crit, U_crit2, U_crit3, U_crit4, U_crit5, U_crit6])
        theta_crit_tot=np.array([thetatUcrit, thetatUcrit2, thetatUcrit3, thetatUcrit4, thetatUcrit5])

        # Database found from internet
        df = pd.read_excel(pathfile+files["database"], usecols=range(11))
        df = df.to_numpy()
        df_name = df[:, -2:-1]
        df = df[:, :-2].astype('float')
        df = np.delete(df, 7, axis=1)
        indAir = [index for index, item in enumerate(df_name.flatten()) if 'WK' in item]
        # Plot
        colors = {
            'Cilinder': 'lightgreen',  # Changed yellow to light green
            'B250': 'blue',
            'Mark1': 'red',
            'Mine GU': 'orange',
            'Mine GY': 'purple',
            'G50': 'cyan',
            'B120': 'magenta',
            'Mark16': 'lime',
            'MK I-IV': 'pink',
            'Mine UMB': 'brown',
            'Mine GC': 'navy',
            'B250 WK': 'olive',
            'Mine GU WK': 'teal',
            'Bol': 'yellow',
            'Sediment': 'black',
            'Cyl': 'purple',
            'Cylin': 'pink',
            'Cilinder2': 'blue',
            'Spher': 'orange',
            'Bolt': 'olive'
        }

        def plot_data_points(df_name, term1, term2):
            """Plot the data points and manage the legend entries."""
            added_entries = set()  # To track added legend entries

            for name in np.unique(df_name):
                indices = np.where(df_name == name)[0]

                # Plot the data points
                plt.loglog(term1[indices], term2[indices], color=colors.get(name, 'black'), marker='*', linestyle='none')  # Data points

                # Add legend entry only if not already added
                if name not in added_entries:
                    plt.plot([], [], color=colors.get(name, 'black'), marker='*', linestyle='none', label=name)
                    added_entries.add(name)

        def model_func(Re, alpha, beta):
            return (alpha*Re ** beta)
        def submersion_depth(d_top, d_bed):
            # Empirical constant for close-packed arrangements
            k = 0.3
            # Calculate depth of submersion
            depth_sub = k * d_top * (d_bed / d_top)
            return depth_sub

        def corrected_velocity(u_top, d_top, d_bed, z_0):
            # Calculate depth of submersion
            depth_sub = submersion_depth(d_top, d_bed)
            # Calculate corrected velocity
            return u_top * (np.log((d_top - depth_sub) / z_0) / np.log(d_top / z_0)), depth_sub
        # Terms
        dbed = df[:, 5] *dbedtok # Rennie
        zbfold = df[:, 6]
        depth_sub = submersion_depth(df[:, 0], dbed)
        zbf = np.maximum(zbfold, depth_sub)
        indx = np.where(zbf!=zbfold)[0]
        u_topold = df[:, 7]
        uchange = (np.log((df[indx, 0] - depth_sub[indx]) / (dbed[indx]/30)) / np.log(df[indx, 0] / (dbed[indx]/30)))
        u_top = u_topold.copy()
        u_top[indx] = uchange * u_top[indx]
        zbsf = zbf/df[:,0]
        D_avgf =  np.sqrt(4 * df[:, 4] / (np.pi * df[:, 1]))

        CM = u_top ** 2 / (g *D_avgf* (df[:, 3] / rho_w - 1))
        fact= np.reshape(1 + 1.5*(df_name=='Cilinder') +0.75*(df_name !='Bol')*(df_name !='Sediment')*(df_name!='Cilinder') , [len(D_avgf),])
        # fact=8.1744*D_avgf/df[:,0] -5.6744
        Ref = (u_top * df[:, 0] * (1-zbsf) /(mu/rho_w))*fact
        # Ref[indAir] = (u_top[indAir] * df[indAir, 0] * (1-zbsf[indAir]) /(12.71e-6))*fact[indAir]
        MFf = rho_w * g * df[:,0]* df[:,4]*(df[:, 3] - rho_w)/(2*mu**2*df[:,1])
        # MFf[indAir] = g * df[indAir,0]* df[indAir,4]*(df[indAir, 3] - 1)/(2*(12.71e-6)**2*df[indAir,1])
        # FFf = np.pi * dbed * rho_w * g * df[:, 0]**2 * D_avgf * (df[:, 3] - rho_w) / (2*mu**2)
        Dk= df[:, 0]/(zbf/0.3)
        FFf = rho_w*df[:, 0]**2/(D_avgf*df[:, 1]*mu**2)*(dbed/df[:, 0])*df[:,4]*(df[:,3]-rho_w)*g

        Re = U_crit_tot *Dmax*(1-zbs)/(mu/rho_w) * dragCoeff  # zbs zbfonly when substantially larger than bed particles: sediment, yes, and Bols maybe
        MF = rho_w*g*Dmax*V_obj*(rho_obj-rho_w)/(2*mu**2*L) * np.ones_like(U_crit_tot)
        FF = rho_w*Dmax**2/(D_avg*L*mu**2)*k*V_obj*(rho_obj-rho_w)*g* np.ones_like(U_crit_tot)

        initial_params = [0.025, 2.01]  #
        lower_bounds = [0, 1.99]
        upper_bounds = [5, 2.01]

        params, covariance = curve_fit(model_func, Ref, MFf, p0=initial_params,  bounds=(lower_bounds, upper_bounds), maxfev=10000 ) # Increase the maximum number of function evaluations)
        beta_fit = 2
        alpha_fit = params[0]
        MF_fit = model_func(Ref, alpha_fit, beta_fit)
        residuals = MFf - MF_fit
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((MFf- np.mean(MFf))**2)
        r_squared = 1 - (ss_res / ss_tot)
        rmse = np.sqrt(np.mean(residuals**2))

        # Plots
        plt.figure(figsize=(10, 6))  # You can adjust the size as needed
        if pp == 1:
            plt.loglog(Dmax/k, Re[0, :], 'rs', label='U_crit (Menzel et al)')  # Red triangles for U_crit
            plt.loglog(Dmax/k, Re[1, :], 'gs', label='U_crit (Rennie anl. et al)')  # Green circles for U_crit2
            plt.loglog(Dmax/k, Re[2, :], 'ms', label='U_crit (Friederichs et al)')  # Magenta squares for U_crit3
            plt.loglog(Dmax/k, Re[3, :], 'ys', label='U_crit (Rennie et al)')  # Yellow crosses for U_crit4
            plt.loglog(Dmax/k, Re[4, :], 'cs', label='U_crit (Chu et al)')  # Yellow crosses for U_crit5
            plt.loglog(Dmax/k, Re[5, :], 's', color='blue',label='U_crit (Deltares)')  # Yellow crosses for U_crit6
        plot_data_points(df_name, Dk, Ref)
        plt.legend(title="Legend", loc='best')  # loc='best' places the legend at the optimal position
        plt.ylabel('Reynolds getal', fontsize=12)  # Font size for better readability
        plt.xlabel('Genormaliseerde Diameter (D/k)', fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(outputpath+"Re_Dk_Testcase_"+str(ii)+".png", dpi=300, bbox_inches='tight')
        # plt.show()
        plt.close()

        plt.figure(figsize=(10, 6))  # You can adjust the size as needed
        if pp == 1:
            indx = (k/Dmax > threshold)
            plt.loglog(Re[0, indx], MF[0, indx], 'rs', label='U_crit (Menzel et al)')  # Red triangles for U_crit
            plt.loglog(Re[1, indx], MF[1, indx], 'gs', label='U_crit (Rennie anl. et al)')  # Green circles for U_crit2
            plt.loglog(Re[2, indx], MF[2, indx], 'ms', label='U_crit (Friederichs et al)')  # Magenta squares for U_crit3
            plt.loglog(Re[3, indx], MF[3, indx], 'ys', label='U_crit (Rennie et al)')  # Yellow crosses for U_crit4
            plt.loglog(Re[4, indx], MF[4, indx], 'cs', label='U_crit (Chu et al)')  # Yellow crosses for U_crit5
            plt.loglog(Re[5, indx], MF[5, indx], 's', color='blue', label='U_crit (Deltares)')  # Yellow crosses for U_crit6

        indx = (1/Dk > threshold)
        plot_data_points(df_name[indx], Ref[indx], MFf[indx])
        MF_upper = MF_fit * (1 + percent_bound)
        MF_lower = MF_fit * (1 - percent_bound)
        plt.loglog(Ref, MF_fit, label=f'Fit: MF = {alpha_fit:.2f} * (Re)^{beta_fit:.2f}', color='blue')
        # plt.loglog(Ref, MF_upper, label='Upper Bound', linestyle='--', color='green')
        # plt.loglog(Ref, MF_lower, label='Lower Bound', linestyle='--', color='orange')
        # plt.fill_between(Ref, MF_lower, MF_upper, color='lightgray', alpha=0.5, label='Bounds Area')
        plt.legend(title="Legend", loc='best')  # loc='best' places the legend at the optimal position
        plt.ylabel('Moment factor', fontsize=12)  # Font size for better readability
        plt.xlabel('Reynolds getal', fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(outputpath+"MF_Re_Testcase_"+str(ii)+".png", dpi=300, bbox_inches='tight')
        # plt.show()
        plt.close()

        plt.figure(figsize=(10, 6))  # You can adjust the size as needed
        if pp == 1:
            indx = (Dmax/k > threshold)
            plt.loglog(Re[0, indx], FF[0, indx], 'rs', label='U_crit (Menzel et al)')  # Red triangles for U_crit
            plt.loglog(Re[1, indx], FF[1, indx], 'gs', label='U_crit (Rennie anl. et al)')  # Green circles for U_crit2
            plt.loglog(Re[2, indx], FF[2, indx], 'ms', label='U_crit (Friederichs et al)')  # Magenta squares for U_crit3
            plt.loglog(Re[3, indx], FF[3, indx], 'ys', label='U_crit (Rennie et al)')  # Yellow crosses for U_crit4
            plt.loglog(Re[4, indx], FF[4, indx], 'cs', label='U_crit (Chu et al)')  # Yellow crosses for U_crit5
            plt.loglog(Re[5, indx], FF[5, indx], 's', color='blue', label='U_crit (Deltares)')  # Yellow crosses for U_crit6
        indx = (Dk > threshold)
        plot_data_points(df_name[indx], Ref[indx], FFf[indx])
        # MF_upper = MF_fit * (1 + percent_bound)
        # MF_lower = MF_fit * (1 - percent_bound)
        # plt.loglog(Ref, MF_fit, label=f'Fit: MF = {alpha_fit:.2f} * (Re)^{beta_fit:.2f}', color='blue')
        # plt.loglog(Ref, MF_upper, label='Upper Bound', linestyle='--', color='green')
        # plt.loglog(Ref, MF_lower, label='Lower Bound', linestyle='--', color='orange')
        # plt.fill_between(Ref, MF_lower, MF_upper, color='lightgray', alpha=0.5, label='Bounds Area')
        plt.legend(title="Legend", loc='best')  # loc='best' places the legend at the optimal position
        plt.ylabel('Kracht factor', fontsize=12)  # Font size for better readability
        plt.xlabel('Reynolds getal', fontsize=12)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(outputpath+"FF_Re_Testcase_"+str(ii)+".png", dpi=300, bbox_inches='tight')
        # plt.show()
        plt.close()

        # Print results
        print("For testcase "+str(ii)+" locatie " + str(header[ii])+":")
        for jj in range(0, len(dt_name)):
            print(dt_name[jj], ':  ', U_crit_tot[-1, jj].T)
        column.append(str(header[ii])+ " [m/s]")
        print(" ")
        minSnel.append(U_crit_tot[-1,:].T)

    minSnel = np.array(minSnel)
    df = pd.DataFrame(np.transpose(minSnel), columns=column)
    df.to_excel(outputxcl, index=False)