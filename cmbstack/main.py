import numpy as np

import maps, stacking


class StackingPipeline:
    """End-to-end stacking pipeline.

    Construct from a power spectrum (:meth:`from_cl`) or from an existing map
    (:meth:`from_map`), then call :meth:`run`.

    Parameters
    ----------
    sky_map : numpy.ndarray
        The HEALPix map to stack on.
    nside : int
        Resolution parameter of the map.

    Attributes
    ----------
    normalized : numpy.ndarray or None
        Set after run(); the normalized map.
    """

    def __init__(self, sky_map, nside):
        self.map = sky_map
        self.nside = nside
        self.normalized = None
        self.peaks = None
        self.patches = None
        self.stacked = None


    @classmethod
    def from_cl(cls, cl_path, nside=128, seed=None):
        """Build a pipeline by simulating a map from a power-spectrum file."""
        # Function 1
        cl = maps.load_cl(cl_path)

        # Function 2
        sky_map = maps.simulate_map(cl, nside, seed)
        return cls(sky_map, nside)
    

    def run(self, size_deg=10.0, reso_arcmin=3.0, profile=True):
        """Run the full stacking loop.

        Parameters
        ----------
        size_deg, reso_arcmin : float
            Patch geometry.
        profile : bool
            Whether to also compute the radial profile.

        Returns
        -------
        result : cmbstack.stack.StackResult
        """

        # Function 3
        self.normalized = maps.normalize_map(self.map)
        
        # Function 4
        self.peaks = stacking.find_hottest_spots(self.normalized)

        # Function 5
        self.patches = stacking.cut_circles_around_spots(self.normalized,self.peaks,self.nside)

        # Function 6
        self.stacked = stacking.average_all_circles(self.patches)

        return self


pipeline = StackingPipeline.from_cl('/Users/isaac/Documents/cmbstack/data/base_plikHM_TTTEEE_lowl_lowE_lensing.minimum.theory_cl', 1024, 42)
pipeline.run()

import matplotlib.pyplot as plt

plt.figure(figsize=(6, 5))
plt.hist(pipeline.stacked, bins=30, edgecolor='black', color='skyblue')
plt.title("The Stacked Peak! (Average temperature of all hot spots)")
plt.xlabel("Temperature value (normalized)")
plt.ylabel("Number of pixels")
plt.grid(True, alpha=0.3)
plt.show()

print("Check the pop-up window for your plot. Success!")