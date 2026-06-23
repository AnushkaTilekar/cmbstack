"""Map normalization: remove the monopole and scale to unit variance."""

import healpy as hp


def normalize_map(m, remove_monopole=True):
    """Subtract the monopole and divide by the standard deviation.

    After this, peak thresholds can be expressed in units of sigma, which is the
    natural convention for peak statistics.

    Parameters
    ----------
    m : array_like
        Input HEALPix map. May contain UNSEEN/NaN pixels, which are ignored in
        the mean and standard deviation.
    remove_monopole : bool, optional
        If True, subtract the mean (monopole) before scaling. Default True.

    Returns
    -------
    m_norm : numpy.ndarray
        The normalized map, with (if monopole removed) mean ~0 and std ~1.
    """
    return hp.remove_monopole(map) / map.std()
