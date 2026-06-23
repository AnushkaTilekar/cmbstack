# FUNCTION 4: Find the hottest spots (peaks) - Detecting the Hottest Local Maxima

#We utilize Healpy's `hotspots` algorithm to identify local maxima. 
#A pixel is considered a local maximum if its temperature value is strictly greater than all its immediate neighboring pixels in the HEALPix grid.

#Once all maxima are found, we rank them by their normalized temperature $T_{\text{norm}}$. We select the top $N=100$ hottest spots:
#\[
#\text{Peak Set} = \{ \hat{n}_p \in \text{Maxima} \mid T_{\text{norm}}(\hat{n}_p) \text{ is in the top } N \}
#\]
#We focus on the hottest peaks because they represent the most significant structures in the Gaussian field and yield the highest signal-to-noise ratio in the final stack.


def find_hottest_spots(cleaned_sky_map, number_of_peaks_to_find=100):
    """
    Finds all the local hot spots, picks the brightest ones, 
    and returns their pixel locations.
    """

    # This finds every single hot spot on the map
    everything, cold_pixels, all_hot_pixels = hp.hotspots(cleaned_sky_map)
    
    # Check the temperature at each hot spot
    temperatures_at_hot_spots = cleaned_sky_map[all_hot_pixels]
    
    # Sort them so the hottest comes first
    order_from_hottest_to_coldest = np.argsort(temperatures_at_hot_spots)[::-1]
    
    # Keep only the top ones (e.g., top 100)
    hottest_peak_pixels = all_hot_pixels[order_from_hottest_to_coldest][:number_of_peaks_to_find]
    
    return hottest_peak_pixels

#-------------------------------------------------------------------------------------
# FUNCTION 5: Cut out little circles around those spots - Extracting Patches (Cutting Circles Around Peaks)

#For each selected peak located at direction vector $\hat{n}_p$, we define a circular aperture (patch) with a fixed angular radius $\Theta = 5^\circ$. 

#To extract the patch, we find every pixel direction $\hat{n}$ on the sphere that falls within this radius. The angular distance $\Delta \theta$ between the peak center and a given pixel is computed via the dot product:
#\[
#\Delta \theta = \arccos\left( \hat{n}_p \cdot \hat{n} \right)
#\]
#We include all pixels satisfying $\Delta \theta \leq \Theta$. This gives us a set of temperature values $T_{\text{norm}}(\hat{n})$ surrounding each peak. Conceptually, we are recentering these patches so that the peak lies exactly at the origin of our local coordinate system.


def cut_circles_around_spots(cleaned_sky_map, hottest_peak_pixels, sky_resolution, circle_size_in_degrees=5):
    """
    For each hot spot, draws a circle (like a cookie cutter) 
    and saves all the temperatures inside that circle.
    """

    # Convert degrees to radians (because Healpy speaks radians)
    circle_size_in_radians = np.radians(circle_size_in_degrees)
    
    all_circles = []  # This will hold all the little cut-out circles
    
    # Loop through each hot spot
    for pixel_location in hottest_peak_pixels:
        # Get the latitude and longitude of this spot
        latitude, longitude = hp.pix2ang(sky_resolution, pixel_location)
        
        # Convert that to a 3D direction arrow (x,y,z)
        center_direction_vector = hp.ang2vec(latitude, longitude)
        
        # Find all the pixel-numbers that fall inside this circle
        pixels_in_circle = hp.query_disc(sky_resolution, center_direction_vector, circle_size_in_radians, inclusive=True)
        
        # Get the actual temperature values at those pixels
        temperatures_in_circle = cleaned_sky_map[pixels_in_circle]
        
        # Save this circle to our list
        all_circles.append(temperatures_in_circle)
    
    return all_circles

#------------------------------------------------------------------------------
# FUNCTION 6: Average all those circles together - STACKING - Averaging the Patches)

#This is the core of our project. To construct the mean stacked peak map, we average the extracted patches pixel-by-pixel. For $N$ peaks, the stacked radial profile is defined as:
#\[
#S(\Delta \theta) = \frac{1}{N} \sum_{i=1}^{N} T_{\text{norm}}^{(i)}(\Delta \theta)
#\]
#where $i$ indexes the selected peaks. 

#By averaging over many independent patches, random fluctuations and noise (which are uncorrelated) average down towards zero, while the coherent, symmetric profile of a typical CMB hot spot remains. This final average gives us the characteristic "peak shape" of the CMB temperature anisotropies.

def average_all_circles(all_circles):
    """
    Takes all the circles, makes them the exact same size, 
    and averages them pixel-by-pixel to create ONE "stacked" circle.
    """

    # Find the smallest circle (just to make sure they all fit together)
    smallest_circle_size = min(len(circle) for circle in all_circles)
    
    # Trim every circle to that exact size
    circles_same_size = [circle[:smallest_circle_size] for circle in all_circles]
    
    # Turn the list into a table (rows = circles, columns = pixels)
    table_of_circles = np.array(circles_same_size)
    
    # Average down the rows (axis=0 means "average each column across all rows")
    average_circle = np.mean(table_of_circles, axis=0)
    
    return average_circle

#-------------------------------------------------------------------
