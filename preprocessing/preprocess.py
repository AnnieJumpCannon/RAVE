
import cPickle as pickle
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from glob import glob
from matplotlib.ticker import MaxNLocator
from random import shuffle

_huge_error = 1e6



def preprocess(spectrum_path, common_wavelengths, continuum_mask, dr4_row,
    neighbour_row=None, neighbour_path=None, output_path=None, fig=None,
    clobber=False):
    """
    Perform preprocessing steps on a RAVE spectrum so that it can be used with
    The Cannon.

    :param spectrum_path:
        The path of the spectrum.
    """
    output_path = output_path or os.path.splitext(spectrum_path)[0] + ".pkl"
    if os.path.exists(output_path) and not clobber:
        return (None, None, None)

    # Interpolate to common wavelength scale.
    data = np.loadtxt(spectrum_path)
    flux = np.interp(common_wavelengths, data[:, 0], data[:, 1],
        left=np.nan, right=np.nan)

    if neighbour_path is not None:
        neighbour_data = np.loadtxt(neighbour_path)
        neighbour_flux = np.interp(common_wavelengths, neighbour_data[:, 0],
            neighbour_data[:, 1], left=np.nan, right=np.nan)


    # Create error array.
    rms = np.nanstd(flux[continuum_mask])
    error = np.ones_like(flux) * rms
    non_finite = ~np.isfinite(flux)
    error[non_finite] = _huge_error

    # Increase variance based on sky lines?
    # TODO: This will depend on how good the sky subtraction was. Since they
    #       are bright it shouldn't be an issue, but we will see.

    # For error region plotting:
    x2 = np.repeat(common_wavelengths, 2)[1:]
    xstep = np.repeat((common_wavelengths[1:] - common_wavelengths[:-1]), 2)
    xstep = np.concatenate(([xstep[0]], xstep, [xstep[-1]]))
    x2 = np.append(x2, x2.max() + xstep[-1])
    x2 -= xstep/2.

    y_lower = np.repeat(flux - error, 2)
    y_upper = np.repeat(flux + error, 2)

    if dr4_row is not None:
        label = "{0}  {1:.0f}/{2:.2f}/{3:+.2f}  (#{4:.0f})".format(
                dr4_row["NAME"], dr4_row["TeffK"], dr4_row["loggK"],
                dr4_row["__M_H_K"], dr4_row["Fiber"])

    else:
        label = ""


    # Make a figure.
    if fig is None:
        fig, ax = plt.subplots(figsize=(14.65, 4.5625))

        ax.plot(common_wavelengths, flux,
            label=label,
            c='k', drawstyle='steps-mid', zorder=10)

        if neighbour is not None and neighbour_path is not None:
            ax.plot(common_wavelengths, neighbour_flux,
                label="{0}  {1:.0f}/{2:.2f}/{3:+.2f}  (#{4:.0f})".format(
                    neighbour_row["NAME"], neighbour_row["TeffK"], 
                    neighbour_row["loggK"], neighbour_row["__M_H_K"],
                    neighbour_row["Fiber"]),
                c='r', drawstyle='steps-mid', zorder=-1)

        else:
            ax.plot([np.nan], [np.nan], c='r', drawstyle='steps-mid', zorder=-1)

        # Show continuum regions (from the mask provided).
        _ = np.where(np.diff(continuum_mask))[0]
        delta = np.diff(common_wavelengths)[0]
        continuum_regions = common_wavelengths[1 + _].reshape(-1, 2)
        for start, end in continuum_regions:
            ax.axvspan(start - 0.5 * delta, end - 0.5 * delta,
                facecolor="#FFCC00", alpha=0.25, zorder=-20, linewidth=0)

        # Styling.
        ax.axhline(1, c="#666666", linestyle=":", zorder=-10)
        ax.xaxis.set_major_locator(MaxNLocator(8))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        
        ax.set_xlim(common_wavelengths[0], common_wavelengths[-1])
        ax.set_ylim(0, 1.1)
        
        ax.set_xlabel(r"$\lambda$ $(\AA{})$")
        ax.set_ylabel(r"$y_{n}$")

        fig.tight_layout()

    else:
        ax = fig.axes[0]

        # Update figure with this data.
        ax.lines[0].set_data([common_wavelengths, flux])
        ax.lines[0].set_label(label)

        if neighbour is not None and neighbour_path is not None:
            ax.lines[1].set_data([common_wavelengths, neighbour_flux])
            ax.lines[1].set_label(
                "{0}  {1:.0f}/{2:.2f}/{3:+.2f}  (#{4:.0f})".format(
                    neighbour_row["NAME"], neighbour_row["TeffK"], 
                    neighbour_row["loggK"], neighbour_row["__M_H_K"],
                    neighbour_row["Fiber"]))

        else:
            ax.lines[1].set_data([[np.nan], [np.nan]])
            ax.lines[1].set_label("")

        # Remove previous error regions.
        ax.collections[-1].remove()

    # Draw error regions.
    ax.fill_between(x2, y_lower, y_upper, where=np.ones_like(x2), 
        color="#AAAAAA", zorder=-1)

    # Update any legend.
    ax.legend(loc="lower right", frameon=False)

    # Create inverse variance arrays.
    ivar = 1.0/(error**2)
    ivar[non_finite] = 0

    # Write the output spectrum to disk.
    with open(output_path, "wb") as fp:
        pickle.dump((flux, ivar), fp, -1)

    return (flux, ivar, fig)


if __name__ == "__main__":

    # Common wavelength scale:
    # 8413.30 to 8790.04 with 0.375 pixels
    # This goes from the 1-percentile of bluest and reddest ends of the spectra.
    wavelengths = np.arange(8413.30, 8790.04 + 0.375, 0.375)

    # Load defined continuum regions.
    continuum_regions_path = "sample3.list"
    with open(continuum_regions_path, "r") as fp:
        regions = [map(float, l.strip().split(":")) for l in fp.readlines()]

    # Build a continuum mask.
    continuum_mask = np.zeros(wavelengths.size, dtype=bool)
    for start, end in regions:
        continuum_mask[(end >= wavelengths) * (wavelengths >= start)] = True

    # Get all the files.
    paths = glob("../reduced-spectra/20??/20*/*.txt")
    shuffle(paths)
    N, fig = len(paths), None

    # Load the last data release.
    dr4 = fits.open("RAVE-DR4.fits")[1].data

    for i, path in enumerate(paths):
        print(i, N, path)

        # Show a similar spectrum.
        name = "{0}_{1}_{2}".format(
            path.split("/")[-2], path.split("/")[-1].split(".")[0],
            path.split(".")[-4])


        index = np.where(dr4["NAME"] == name)[0]

        if len(index) > 0:

            dr4_row = dr4[index[0]]

            if np.isfinite(dr4["TeffK"][index]):
                distances = (((dr4["TeffK"] - dr4["TeffK"][index])/1000)**2 \
                            + (dr4["loggK"] - dr4["loggK"][index])**2 \
                            + (dr4["__M_H_K"] - dr4["__M_H_K"][index])**2)**0.5

                # Check this star has parameters?
                distances[~np.isfinite(distances)] = 10e6

                closest_indices = np.argsort(distances)

                for closest_index in closest_indices[1:]:
                    neighbour_path = glob(
                        "../reduced-spectra/20??/{0}/{1}.*.{2}.*.txt".format(
                            *tuple(dr4["NAME"][closest_index].split("_"))))

                    if len(neighbour_path) > 0 \
                    and os.path.exists(neighbour_path[0]): break

                neighbour_path = neighbour_path[0]

                neighbour = dr4[closest_index]

            else:
                neighbour, neighbour_path = None, None

        else:
            dr4_row, neighbour, neighbour_path = None, None, None

        # Do the pre-processing.
        flux, ivar, fig = preprocess(path, wavelengths, continuum_mask,
            dr4_row=dr4_row,
            neighbour_row=neighbour, neighbour_path=neighbour_path,
            fig=fig)
        if fig is None: continue

        # Save the figure.
        figure_path = "../quicklook/{0}.png".format(
            os.path.splitext(path)[0][len("../reduced-spectra/"):])
        if not os.path.exists(os.path.dirname(figure_path)):
            os.makedirs(os.path.dirname(figure_path))

        fig.savefig(figure_path, dpi=150)


