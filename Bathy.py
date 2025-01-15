
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from sklearn.decomposition import PCA
from scipy.interpolate import interp1d, NearestNDInterpolator
from scipy.spatial import cKDTree

from scipy.optimize import curve_fit

def cos_func(x, A, omega):
    return A * np.cos(omega * x)

# Load the data with pandas, specifying comma as the delimiter
def bathy(files):
    # Input
    pathfile = ".\\Excel\\"
    outputpath = ".\\Plaatjes\\"
    outputxcl = pathfile + files["output_file"]

    # Testcases
    dk = pd.read_excel(pathfile + files["testcases"], header=None)
    dk = dk.to_numpy()
    columns = dk[1:, 1]
    header = dk[1:, -1]
    indx = np.unique(header, return_index=True)
    indx = indx[1]
    header = header[indx].astype("str")
    columns = columns[indx].astype("str")
    columns = np.core.defchararray.add(columns, " [m]")

    # Settings
    maxsiz = 1600000
    num_levels = 50
    width = 30

    # OOlijst
    dt = pd.read_excel(pathfile+files["OOlijst"])
    dt = dt.to_numpy()
    dt_name = dt[:, -3:-2]
    dt = dt[:, :-3].astype('float')
    lenvect = dt[:, 1]
    depthTotTot = np.reshape(dt_name, [12,])
    column = np.hstack([['Object'], columns])

    for jj in range(0, len(header)):
        file_path = ".\\DataBodem\\" + header[jj]
        print(file_path)
        data_loaded = np.load(file_path)
        X = data_loaded['X']
        Y = data_loaded['Y']
        Z = data_loaded['Z']
        fact  = round(8000000*0.001)
        X = X[::fact]
        Y = Y[::fact]
        Z = Z[::fact] - np.mean(Z)

        # Fit a spline to the X, Y coordinates (parametric fit)
        X = data_loaded['X']
        Y = data_loaded['Y']
        Z = data_loaded['Z'] - np.mean(Z)

        coords = np.column_stack((X, Y))
        pca = PCA(n_components=2)
        coords_pca = pca.fit_transform(coords)
        indx = np.argsort(coords_pca[:, 0])
        main_axis = coords_pca[indx, 0]  # X-like coordinates in PCA space
        perpendicular_axis = coords_pca[indx, 1]  # Y-like coordinates in PCA space
        Z = Z[indx]
        degree = 3
        coeffs = np.polyfit(main_axis, perpendicular_axis, degree)
        vect =  np.linspace(min(main_axis), max(main_axis), round(abs(np.max(main_axis) - np.min(main_axis))/0.01))
        poly_curve = np.polyval(coeffs, vect)
        centerline_points = np.column_stack((vect, poly_curve))
        centerline_tree = cKDTree(centerline_points)
        points = np.column_stack([main_axis, perpendicular_axis])
        distances, _ = centerline_tree.query(points)
        channel_width_half = width/2  # Half of 100 meters
        channel_indices = np.where(distances <= 0.025)[0]
        interpolator = NearestNDInterpolator(list(zip(main_axis[channel_indices], perpendicular_axis[channel_indices])), Z[channel_indices])
        Xn = centerline_points[:, 0]
        Yn  = centerline_points[:, 1]
        Zn = interpolator(Xn, Yn)

        indx = np.argsort(Xn)
        Xn = Xn[indx]
        Yn = Yn[indx]
        Zn= Zn[indx]
        channel_indices = np.where(distances <= channel_width_half)[0]
        X = main_axis[channel_indices]
        Y = perpendicular_axis[channel_indices]
        Z = Z[channel_indices]
        Z = Z - np.mean(Z)
        X = X-np.mean(X)
        Y = Y-np.mean(Y)
        distance = np.sqrt(X**2 + Y**2)
        indx = np.argsort(distance)
        X = X[indx[:maxsiz]]
        Y = Y[indx[:maxsiz]]
        Z = Z[indx[:maxsiz]]

        # Save figure
        plt.figure()
        levels = np.linspace(np.min(Z), np.max(Z), num_levels)
        contour = plt.tricontourf(X, Y, Z, levels=levels, cmap='jet')

        plt.plot(centerline_points[:, 0], centerline_points[:, 1], 'b--')
        plt.colorbar(contour)
        plt.title("Map")
        plt.xlabel("X in m")
        plt.ylabel("Y in m")
        plt.xlim([np.min(X), np.max(X)])
        plt.ylim([np.min(Y), np.max(Y)])
        plt.savefig(file_path + ".png", dpi=300, bbox_inches='tight')

        # Fourier
        dx = np.diff(Xn)  # Difference in X direction
        dy = np.diff(Yn)  # Difference in Y direction
        ds = np.sqrt(dx**2 + dy**2)
        ds = np.cumsum(ds)
        ds = np.insert(ds, 0, 0)
        ds, indices = np.unique(ds, return_index=True)
        Zn = Zn[indices]
        Yn = Yn[indices]
        Xn = Xn[indices]
        ds_uniform = np.linspace(np.min(ds), np.max(ds), round(abs(np.max(ds) - np.min(ds))/0.01))
        interpolator = interp1d(ds, Zn, kind='cubic', fill_value="extrapolate")
        Zn = interpolator(ds_uniform)
        ds = ds_uniform
        Znorg = Zn - np.mean(Zn)
        Zn_fft = np.fft.fft(Znorg)

        sampling_interval = ds[1] - ds[0]  # assuming ds is the arc length increment
        frequencies = np.fft.fftfreq(len(Zn), d=sampling_interval)

        positive_frequencies = frequencies[:len(frequencies) // 2]
        positive_Zn_fft = Zn_fft[:len(Zn_fft) // 2]

        magnitude = np.abs(positive_Zn_fft)
        depthTot = []
    #---
        for ii in range(0, len(lenvect)):
            cutoff_frequency = 1/(2 * lenvect[ii])* 0.5
            cutoff_frequency2 = 1/(0.5 * lenvect[ii]) * 0.5
            # cutoff_frequency2 = 1/(1e-6)* 0.5
            lowpass_filter = np.abs(frequencies) > cutoff_frequency
            lowpass_filter2 = np.abs(frequencies) < cutoff_frequency2

            Zn_fft_filtered = Zn_fft * lowpass_filter * lowpass_filter2
            Zn_filtered = np.fft.ifft(Zn_fft_filtered)
            Zn_filtered = np.real(Zn_filtered)

            # Find zero-crossings
            zero_crossings = np.where(np.diff(np.sign(Zn_filtered)))[0]  # Find zero-crossings
            upcrossings = zero_crossings[Zn_filtered[zero_crossings] < 0]

            # Initialize a list to store the depths (max peak-to-peak difference)
            depth = []

            # Iterate over each pair of consecutive zero-crossings
            nn = 0
            val = 3
            for i in range(val, len(upcrossings) - 1-val):
                # Define the segment between two consecutive zero-crossings (full wave)
                start = upcrossings[i]
                end = upcrossings[i + 1]+2
                # period = ds[end]-ds[start]
                # omega = 2 * np.pi / period
                # popt, _ = curve_fit(lambda x, A: cos_func(x, A, omega),  ds[start:end], Zn_filtered[start:end], p0=[ np.max(Zn_filtered[start:end])], maxfev=20000)

                # Find the minimum value in this segment
                segment_min = abs(np.min(Zn_filtered[start:end]))
                segment_max = abs(np.max(Zn_filtered[start:end]))

                if i > len(upcrossings)/2 and nn==0:
                    nn = 1
                    named = np.char.replace(header[jj], ".npz", "")
                    plt.figure(figsize=(10, 6))
                    plt.plot(ds[start:end], Zn_filtered[start:end], color='red')
                    plt.xlabel('Afstand (m)')
                    plt.ylabel('Diepte gefilterd voor lengte (m)')
                    plt.title("Testcase_" + str(named) + "_Objectlength_" + str(lenvect[ii]) )
                    plt.tight_layout()
                    plt.savefig(outputpath + "Testcase_" + str(named) + "_Objectlength_" + str(lenvect[ii]) + "_segment.png", dpi=300, bbox_inches='tight')
                    plt.close()

                # Calculate the depth as the difference between max and min
                depth.append((segment_min+segment_max)/2)

            # depth = np.mean(Zn_filtered[indx1]) - np.mean(Zn_filtered[indx2])
            depth = np.mean(np.array(depth))
            depthTot.append(depth)
            print('Diepte: ', depth, ' met lengte ', lenvect[ii])
    #----
        depthTotTot = np.vstack([depthTotTot, np.array(depthTot)])
        cutoff_frequency = 1 / (1.1 * max(lenvect)) * 0.5
        cutoff_frequency2 = 1 / (0.9 * min(lenvect)) * 0.5
        # cutoff_frequency2 = 0
        lowpass_filter = np.abs(frequencies) > cutoff_frequency
        lowpass_filter2 = np.abs(frequencies) < cutoff_frequency2
        Zn_fft_filtered = Zn_fft * lowpass_filter * lowpass_filter2
        Zn_filtered = np.fft.ifft(Zn_fft_filtered)
        Zn_filtered = np.real(Zn_filtered)

        indx1, _ = find_peaks(Zn_filtered)
        indx2, _ = find_peaks(-(Zn_filtered))

        depth = np.mean(Zn_filtered[indx1]) - np.mean(Zn_filtered[indx2])

        print('Depth: ', depth, ' total average')
        Zn_filtered_fft = np.abs(np.fft.fft(Zn_filtered)[:len(np.fft.fft(Zn_filtered)) // 2])

        # Plot the original data
        plt.figure(figsize=(10, 6))
        plt.subplot(2, 1, 1)
        plt.plot(ds,  Znorg, label="Originele Data")
        plt.plot(ds,  Znorg - Zn_filtered ,color='red', label="Low pass Data")
        plt.plot(ds,  Zn_filtered ,color='orange', label="Low pass Data")
        plt.xlabel("Afstand in m")
        plt.ylabel("Relatieve hoogte bodem in m")
        plt.title("Originele bodem data")
        plt.grid(True)

        # Plot the magnitude of the Fourier transform (frequency domain)
        plt.subplot(2, 1, 2)
        plt.plot(positive_frequencies, magnitude, label="Fourier Transform Waardee", color='blue')
        plt.plot(positive_frequencies, Zn_filtered_fft,  label="Fourier Transform Waarde", color='orange')
        plt.xlabel("Frequentie (1/m)")
        plt.ylabel("Waarde")
        plt.title("Fourier Transform van bodem hoogte")
        plt.grid(True)
        # plt.xlim([0,5])

        # Save the figure
        plt.tight_layout()
        named = np.char.replace(header[jj], ".npz", "")
        plt.savefig(outputpath +str(named)+ "_fourier_transform.png", dpi=300, bbox_inches='tight')
        plt.close()

    df = pd.DataFrame(np.transpose(depthTotTot), columns=column)
    df.to_excel(outputxcl, index=False)