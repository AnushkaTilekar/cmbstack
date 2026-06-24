# CMB Peak Stacking Pipeline

## Overview

This repository contains our Python implementation for stacking patches of the Cosmic Microwave Background (CMB) temperature sky. Our goal is to take a theoretical power spectrum, generate a synthetic sky map, locate the hottest local maxima (peaks), and average the surrounding patches. This averaging process, known as "stacking," enhances the symmetric signal of the typical peak profile while suppressing random noise and cosmic variance.

Below, we document the mathematical reasoning behind each of our six core functions.

## Dependencies

- `healpy` — handles the round sky maps
- `numpy` — handles numbers and math
- `matplotlib` — draws pictures

## Summary of Workflow

1. **Input:** Theoretical $D_\ell$ spectrum.
2. **Convert:** Change $D_\ell$ to $C_\ell$.
3. **Simulate:** Generate a random Gaussian CMB sky map.
4. **Standardize:** Remove the mean and normalize by $\sigma$.
5. **Identify:** Find and rank the local maxima.
6. **Extract:** Cut out $5^\circ$ discs around the top 100 peaks.
7. **Average:** Stack the patches to reveal the universal peak profile.

---

## Functions

### 1. Reading the Power Spectrum (Converting $D_\ell$ to $C_\ell$)

**Function:** `read_power_spectrum_from_file`

The data provided by our supervisors contains the power spectrum in the form of $D_\ell^{TT}$, defined as:

$$D_\ell \equiv \frac{\ell(\ell+1)}{2\pi} C_\ell$$

However, the Healpy library expects the angular power spectrum $C_\ell$ as input for generating maps. Therefore, we rearrange the formula to convert our data:

$$C_\ell = \frac{2\pi}{\ell(\ell+1)} D_\ell$$

We set $C_0 = 0$ to avoid the division-by-zero singularity at $\ell=0$, as the monopole contributes no fluctuating signal.

**Implementation notes:**
- Opens the file, skipping the first line (which has the `#` and headers).
- Column 0: the scale of the bumps.
- Column 1: the Temperature power values (TT) — the only column we use.
- Converts using: `Cl = Dl * 2 * pi / (l * (l + 1))`.
- Sets the first value to 0 to avoid division-by-zero.

---

### 2. Generating the Synthetic CMB Map

**Function:** `create_fake_sky_map`

We generate a random, Gaussian-distributed realization of the CMB sky. Healpy achieves this by drawing spherical harmonic coefficients $a_{\ell m}$ from a Gaussian distribution with variance $C_\ell$, formally defined as:

$$\langle |a_{\ell m}|^2 \rangle = C_\ell$$

The temperature at any point on the sphere $\hat{n}$ is then computed via the spherical harmonic transform:

$$T(\hat{n}) = \sum_{\ell=0}^{\ell_{\text{max}}} \sum_{m=-\ell}^{\ell} a_{\ell m} \, Y_{\ell m}(\hat{n})$$

where $Y_{\ell m}$ are the spherical harmonic functions. The `synfast` routine performs this exact computation to produce our mock sky map. `sky_resolution=128` gives a decent quality picture.

---

### 3. Preprocessing the Map (Removing Monopole & Normalizing)

**Function:** `remove_average_and_scale_map`

Before searching for peaks, we must standardize our map to ensure the detection threshold is uniform. First, we remove the monopole (the average temperature) from the map:

$$T_{\text{fluctuation}}(\hat{n}) = T(\hat{n}) - \langle T \rangle$$

Next, we divide by the standard deviation $\sigma$ to make the map dimensionless and unit-variant:

$$T_{\text{norm}}(\hat{n}) = \frac{T_{\text{fluctuation}}(\hat{n})}{\sigma}$$

where $\sigma = \sqrt{\langle (T - \langle T \rangle)^2 \rangle}$. After this step, our map has a mean of zero and a standard deviation of one, meaning peaks are measured in units of "sigma" above the average.

---

### 4. Detecting the Hottest Local Maxima

**Function:** `find_hottest_spots`

We utilize Healpy's `hotspots` algorithm to identify local maxima. A pixel is considered a local maximum if its temperature value is strictly greater than all its immediate neighboring pixels in the HEALPix grid.

Once all maxima are found, we rank them by their normalized temperature $T_{\text{norm}}$. We select the top $N=100$ hottest spots:

$$\text{Peak Set} = \{ \hat{n}_p \in \text{Maxima} \mid T_{\text{norm}}(\hat{n}_p) \text{ is in the top } N \}$$

We focus on the hottest peaks because they represent the most significant structures in the Gaussian field and yield the highest signal-to-noise ratio in the final stack.

**Implementation notes:**
- Finds every local hot spot on the map.
- Checks the temperature at each hot spot.
- Sorts from hottest to coldest.
- Keeps only the top N (default: 100).

---

### 5. Extracting Patches (Cutting Circles Around Peaks)

**Function:** `cut_circles_around_spots`

For each selected peak located at direction vector $\hat{n}_p$, we define a circular aperture (patch) with a fixed angular radius $\Theta = 5^\circ$.

To extract the patch, we find every pixel direction $\hat{n}$ on the sphere that falls within this radius. The angular distance $\Delta \theta$ between the peak center and a given pixel is computed via the dot product:

$$\Delta \theta = \arccos\left( \hat{n}_p \cdot \hat{n} \right)$$

We include all pixels satisfying $\Delta \theta \leq \Theta$. This gives us a set of temperature values $T_{\text{norm}}(\hat{n})$ surrounding each peak. Conceptually, we are recentering these patches so that the peak lies exactly at the origin of our local coordinate system.

**Implementation notes:**
- Converts degrees to radians (Healpy uses radians internally).
- For each hot spot: gets the latitude/longitude, converts to a 3D direction vector (x, y, z), finds all pixels inside the disc, and saves their temperatures.

---

### 6. Stacking (Averaging the Patches)

**Function:** `average_all_circles`

This is the core of our project. To construct the mean stacked peak map, we average the extracted patches pixel-by-pixel. For $N$ peaks, the stacked radial profile is defined as:

$$S(\Delta \theta) = \frac{1}{N} \sum_{i=1}^{N} T_{\text{norm}}^{(i)}(\Delta \theta)$$

where $i$ indexes the selected peaks.

By averaging over many independent patches, random fluctuations and noise (which are uncorrelated) average down towards zero, while the coherent, symmetric profile of a typical CMB hot spot remains. This final average gives us the characteristic "peak shape" of the CMB temperature anisotropies.

**Implementation notes:**
- Finds the smallest circle size to ensure all patches align correctly.
- Trims every circle to that exact size.
- Arranges patches into a 2D array (rows = circles, columns = pixels).
- Averages down the columns (`axis=0`) to produce the final stacked result.
