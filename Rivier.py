import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import interp2d, interp1d
import pandas as pd
# Code to predict flow pattern cross sectional which are relatively wide
def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def rivier(files):
    pathfile = ".\\Excel\\"
    outputpath = ".\\Plaatjes\\"
    outputxcl = pathfile+files["output_file"]

    dataTot = pd.read_excel(pathfile+files["testcases"])
    dataOO = pd.read_excel(pathfile+files["snelheid"])
    dataOV = pd.read_excel(pathfile+files["OOlijst"])
    header = dataOO.columns.values[1:]
    dataOO = dataOO.to_numpy()
    dataTot = dataTot.to_numpy()
    dataOV = dataOV.to_numpy()
    no = np.size(dataTot, 0)
    objects = dataOO[:, 0]


    for ii in range(0, no):
        # filename
        fileName = dataTot[ii, 5] # downloaded from bathymetry Nederland for a river
        title = dataTot[ii, 1] + " - " + dataTot[ii, 2]
        hbed = dataTot[ii, 4] / 100 # waterhoogte t.o.v. NAP in m
        Qcalib = dataTot[ii, 3] # Calibration value for the debit
        hu = 0.0 # Depth where maximum velocity is found with positive value

        # Find indices in header_array that contain the search string
        indices = [idx for idx, col in enumerate(header) if dataTot[ii, 1] in col]
        threshold = dataOO[:, indices[0]+1] # movement of transport m/s

        # Make river shape
        data = np.genfromtxt(".\\Dwarsdoorsnede\\"+fileName, delimiter=';', dtype=str, skip_header=1)

        # Initialize a list to store indices of columns containing at least one float
        columns_with_floats = []

        # Iterate over columns starting from column 3 (index 2)
        for col_idx in range(2, data.shape[1]):
            try:
                # Try to convert the column to float
                col_data = data[:, col_idx]
                float_values = [float(val) for val in col_data if val.replace('.', '', 1).isdigit()]

                # If at least one value is a float, store the column index
                if float_values:
                    columns_with_floats.append(col_idx)
            except ValueError:
                # Ignore columns that cannot be converted to floats
                continue
        # Print the indices of columns containing at least one float
        data = data[:, columns_with_floats]

        # Filter rows where all values are valid floats
        data = np.array([row for row in data if all(is_float(value) for value in row)])
        data = data.astype("float")
        data[:, 0] = data[:, 0] - data[0, 0]
        h = data[:, 1] - hbed
        X = data[:,0]

        # Channel abnks
        try:
            indx1 = np.where((h[:-1] > 0) & (h[1:] < 0))[0][0]
        except:
            indx1=0
        try:
            indx2 = np.where((h[:-1] < 0) & (h[1:] > 0))[0][0] +1
        except:
            indx2 = len(h)

        h = h[indx1:indx2+1]
        X = X[indx1:indx2+1]
        hnew=h.copy()
        hnew[[0, -1]] = 0
        X[0] = X[0] + (X[1] - X[0]) * (0 - h[0]) / (h[1] - h[0])
        X[-1] = X[-1 - 1] + (X[-1] - X[-1 - 1]) * (0 - h[-1 - 1]) / (h[-1] - h[-1 - 1])
        h=hnew
        X=X-X[0]
        area = np.trapezoid(abs(h), X)
        xm =  np.trapezoid(abs(h) * X, X) / area
        Um = Qcalib/area

        # indxloc = np.argmin(np.abs(X-(X[-1]- X[0]) * xperc)) # Assumption of maximum velocity lcoation
        indxloc = np.argmin(np.abs(X-xm)) # Assumption of maximum velocity lcoation
        xloc = X[indxloc]
        def equation(M, Phi):
            return (np.exp(M) / (np.exp(M) - 1)) - (1 / M) - Phi
        # Parameters
        M = 0.7
        coeff = 0.55
        tolerance_outer = 1e-3
        tolerance_inner = 1e-4
        max_iterations_outer = 500
        max_iterations_inner = 200
        initial_step_size = 0.01
        min_step_size = 1e-6
        u = 0
        iteration_outer = 0
        umax=0
        # Start the outer convergence loop
        while abs(np.mean(u) - Um) > tolerance_outer and iteration_outer < max_iterations_outer:
            iteration_outer += 1
            residual = abs(np.mean(u) - Um)
            print(f"Iteratie {iteration_outer} | Residu: {residual} | Coeff: {coeff} | Max: {umax}")

            # Adjust `coeff` adaptively based on the residual
            step_size = initial_step_size * residual
            step_size = max(step_size, min_step_size)  # Ensure it doesn't get too small

            # Update `coeff` in the right direction
            if np.mean(u) > Um:
                coeff += step_size  # Increase `coeff` if `u` is too low
            else:
                coeff -= step_size  # Decrease `coeff` if `u` is too high

            umaxsurf = Um / coeff  # Based on Moramarco & Singh 2010
            umax = umaxsurf
            umaxold = 0
            iteration_inner = 0

            # Inner loop for calculating `umax`
            while abs(umax - umaxold) > tolerance_inner and iteration_inner < max_iterations_inner:
                iteration_inner += 1
                umaxold = umax
                delta = abs(h[indxloc]) / (abs(h[indxloc]) - hu)
                umax = umaxsurf / (1 / M * np.log(1 + (np.exp(M) - 1) * delta * np.exp(1 - delta)))
                PhiM = Um / umax
                M = fsolve(equation, M, args=(PhiM))[0]

            # Final adjustment to M
            M = fsolve(equation, M, args=(0.55))[0]

            # Grid setup
            num_grid_points = 100
            X_fine = np.linspace(X.min(), X.max(), num_grid_points)
            h_interpolated = np.interp(X_fine, X, h)
            heights = np.linspace(h_interpolated.min(), h_interpolated.max(), num_grid_points)
            X_grid, H_grid = np.meshgrid(X_fine, heights)
            H_grid = -np.clip(-H_grid, 0, -h_interpolated[:, np.newaxis].T)
            X_locations = X_grid.flatten()
            H_locations = H_grid.flatten()
            locations = np.column_stack((X_locations, H_locations))
            locations = np.unique(locations, axis=0)
            h_interpolated = np.interp(locations[:, 0], X_fine, h_interpolated)

            # Hydraulics calculations
            perimeter = np.sum(np.sqrt(np.diff(X)**2 + np.diff(h)**2))
            Rh = area / perimeter

            # Umax profile for width calculation
            xv = (locations[:, 0] - xloc)
            y = abs(locations[:, 1])
            D = abs(h_interpolated)
            xs = np.where(xv < 0, abs(np.min(xv)), abs(np.max(xv)))
            umaxprof = umax * (1 - (xv / xs)**2) # small width rivers Discussable with Peter
            # umaxprof = umax * (1 - (xv / xs)**2)**0.5 # large width rivers
            u = np.zeros(len(locations))

            # Compute velocities
            for jj in range(len(u)):
                if D[jj] == 0 and hu == 0:
                    u[jj] = 0
                elif y[jj] < D[jj]:
                    delta_y = (D[jj] - y[jj]) / (D[jj] - hu)
                    u[jj] = umaxprof[jj] / M * np.log(1 + (np.exp(M) - 1) * delta_y * np.exp(1 - delta_y))

        # Plot
        print('De ratio diepte vs ruwheid: ', str(np.exp((PhiM -0.51)/0.11)if PhiM < 0.66 else np.inf))
        plt.figure(figsize=(10, 6))
        scatter = plt.scatter(locations[:, 0], locations[:, 1], c=u, marker='o', cmap='viridis')
        plt.colorbar(scatter, label='U (m/s)')  # Colorbar for velocity
        plt.plot(X, h, color='red')
        plt.xlabel('Coordinaat in breedterichting rivier (m)')
        plt.ylabel('Diepte (m)')
        plt.title('Snelheidsverdeling in open kanaal stroming voor dwarsdoorsnede')

        # Save the figure
        plt.tight_layout()
        plt.savefig(outputpath + "Testcase_" + str(dataTot[ii, 1]) + "_Qcalib_" + str(Qcalib) + "_rivier.png", dpi=300,
                    bbox_inches='tight')

        plt.close()

        # at height of object
        for jj in range(0, len(objects)):
            # Specify the desired height D
            u_interpolated = np.zeros_like(np.unique(locations[:, 0]))
            DD = dataOV[jj, 0]
            for kk in range(0, len(np.unique(locations[:, 0]))):
                x_val = np.unique(locations[:, 0])[kk]
                indice = np.where(locations[:, 0] == x_val)[0][0]
                if (D[indice] -DD<= 0):
                    u_interpolated[kk] = 0
                else:
                    delta_y = (DD) / (D[indice] - hu)
                    u_interpolated[kk] = umaxprof[indice] / M * np.log(1 + (np.exp(M) - 1) * delta_y * np.exp(1 - delta_y))

            x_line = np.linspace(min(locations[:, 0]), max(locations[:, 0]), 500)

            plt.figure(figsize=(10, 6))
            plt.plot(np.unique(locations[:, 0]), u_interpolated, color='red')
            plt.plot(x_line, np.ones_like(x_line)*threshold[jj], 'b--',)
            plt.xlabel('Coordinaat in breedterichting rivier (m)')
            plt.ylabel('Snelheid (m/s)')
            plt.title("Testcase_" + str(dataTot[ii, 1]) + "_Qcalib_" + str(Qcalib) + "_" +str(objects[jj]))
            plt.tight_layout()
            plt.savefig(outputpath + "Testcase_" + str(dataTot[ii, 1]) + "_Qcalib_" + str(Qcalib) + "_" +str(objects[jj]) + "_rivier.png", dpi=300,
                    bbox_inches='tight')
            plt.close()
