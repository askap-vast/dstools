def gaussian_fit(df, threshold=1.5, binsize=0, max_iterations=100, min_improvement=0.0005, 
                       window_size=0.1, noise_floor=0.15, verbose=True, errors=True, show_peaks=False):
    """
    Enhanced iterative Gaussian fit for light curves with optimised peak detection and fitting.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input DataFrame containing time and flux columns ('t_mjd' and 'y_i').
    threshold : float, optional
        Threshold multiplier for the standard deviation to consider flux above the baseline. Default is 1.5.
    binsize : int, optional
        Bin size for rebinning the data. Default is 0 (no rebinning).
    max_iterations : int, optional
        Maximum number of peaks to fit before stopping. Default is 100.
    min_improvement : float, optional
        Minimum fractional improvement in residual variance required to continue. Default is 0.0005 (0.05%).
    window_size : float, optional
        Size of the window around each peak in days. Default is 0.1 (~144 minutes).
    noise_floor : float, optional
        Estimated noise floor level - once residuals reach this level, fitting stops. Default is 0.15.
    verbose : bool, optional
        Whether to print progress information. Default is True.

    Returns:
    --------
    results : list of dict
        Each dict contains the parameters of the fitted Gaussian (amplitude, center, width),
        as well as the individual peak model and the cumulative model at each step.
    """
    import numpy as np
    from scipy.signal import find_peaks, savgol_filter
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import time
    import warnings
    
    # Suppress specific RuntimeWarnings from curve_fit that can clutter output
    warnings.filterwarnings("ignore", category=RuntimeWarning, 
                          message=".*covariance of the parameters.*")

    def gaussian(x, amp, center, width):
        """Gaussian function with positive amplitude."""
        return amp * np.exp(-((x - center) ** 2) / (2 * width ** 2))

    def multi_gaussian(x, *params):
        """Multiple Gaussian model for refining fits."""
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            amp = params[i]
            center = params[i+1]
            width = params[i+2]
            y += gaussian(x, amp, center, width)
        return y

    def rebin(array, binsize):
        """Rebin an array to reduce noise."""
        left = len(array) % binsize
        if left != 0:
            array = array[:-int(left)]
        return np.mean(array.reshape(-1, binsize), axis=1)
    
    def estimate_baseline(y, percentile=30):
        """Estimate baseline using a low percentile."""
        return np.percentile(y, percentile)

    # Extract time and flux
    t = df['t_mjd'].values
    y = df['y_i'].values
    e_y = df['dy_i'].values
    if binsize > 0:
        if verbose:
            print(f"Rebinning data with bin size {binsize}")
        t = rebin(t, binsize)
        y = rebin(y, binsize)

    # Apply light smoothing to make peak detection more robust
    try:
        # Use higher window length for better smoothing
        window_len = min(11, len(y) // 5)
        if window_len % 2 == 0:  # Must be odd
            window_len += 1
        y_smooth = savgol_filter(y, window_length=window_len, polyorder=2)
    except:
        # If savgol fails (e.g., data too short), use original data
        y_smooth = y.copy()
        if verbose:
            print("Savgol filter failed, using original data")

    # Estimate and remove baseline for better peak detection
    baseline = estimate_baseline(y_smooth)
    if verbose:
        print(f"Estimated baseline: {baseline:.4f}")
    
    # Working copies
    residuals = y.copy()
    overall_model = np.zeros_like(y)
    results = []
    peak_models = []  # Store models separately
    peak_params = []  # Store all peak parameters for later refinement
    
    # Calculate initial variance for convergence check
    initial_variance = np.var(residuals)
    prev_variance = initial_variance
    
    # Calculate standard deviation for amplitude filtering
    std_dev = np.std(residuals)
    
    if verbose:
        print(f"Starting iterative fitting. Initial variance: {initial_variance:.6f}")
        print(f"Standard deviation: {std_dev:.6f}")
    
    # Create plots only once, will update later
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]}, figsize=(12, 8), sharex=True)
    ax1.plot(t, y, 'k.', label='Raw Data', alpha=0.5)
    
    start_time = time.time()
    
    # Main loop with explicit iteration count
    for iteration in range(max_iterations):
        # Find peaks in the residuals with more sensitive parameters
        current_std_dev = np.std(residuals)
        
        # Use height relative to baseline rather than absolute zero
        peaks, peak_properties = find_peaks(
            residuals, 
            height=threshold * current_std_dev,
            distance=int(0.005 / np.median(np.diff(t))) or 3,  # Dynamic distance based on time sampling
            prominence=0.5 * current_std_dev  # Add prominence requirement
        )
        
        if len(peaks) == 0:
            if verbose:
                print(f"No more peaks found above threshold ({threshold} * {current_std_dev:.6f}). Stopping.")
            break

        # Identify the highest peak
        peak_idx = peaks[np.argmax(residuals[peaks])]
        peak_time = t[peak_idx]
        peak_flux = residuals[peak_idx]
        
        if verbose:
            print(f"Iteration {iteration+1}: Fitting peak at time {peak_time:.6f} with flux {peak_flux:.6f}")
        
        # Select a wider window around the peak for fitting
        window_mask = np.abs(t - peak_time) < window_size
        t_window = t[window_mask]
        y_window = residuals[window_mask]

        if len(t_window) < 7:  # Need more points for a robust fit
            if verbose:
                print(f"Not enough points around peak at {peak_time}. Skipping.")
            # Mark this region to avoid repeated attempts
            residuals[peak_idx] = baseline
            continue

        # Try different initial width guesses and fitting approaches
        width_guesses = [0.005, 0.01, 0.02, 0.001, 0.05]
        success = False
        best_fit = None
        best_rsq = -np.inf
        
        for width_guess in width_guesses:
            # Initial guesses for Gaussian parameters
            amp_guess = peak_flux
            center_guess = peak_time
            
            try:
                # Use broader bounds for better fitting
                popt, pcov = curve_fit(
                    gaussian, 
                    t_window, 
                    y_window, 
                    p0=[amp_guess, center_guess, width_guess],
                    bounds=([0, t_window.min(), 0.0001], 
                            [np.inf, t_window.max(), 0.2]),
                    maxfev=10000,  # Increase max function evaluations
                    method='trf'   # More robust fitting method
                )
                
                amp, center, width = popt
                
                # Calculate R-squared to judge fit quality
                y_pred = gaussian(t_window, amp, center, width)
                ss_res = np.sum((y_window - y_pred) ** 2)
                ss_tot = np.sum((y_window - np.mean(y_window)) ** 2)
                r_squared = 1 - (ss_res / ss_tot)
                
                # Check if the fit makes sense and is better than previous attempts
                if width > 0.0001 and width < 0.2 and amp > 0 and r_squared > best_rsq:
                    best_fit = (amp, center, width)
                    best_rsq = r_squared
                    success = True
                
            except Exception as e:
                if verbose:
                    print(f"  Fit attempt failed with width {width_guess}: {e}")
                continue
        
        if not success:
            if verbose:
                print(f"Failed to fit peak at {peak_time} after trying multiple width guesses.")
            # Mark this region to avoid repeated attempts
            residuals[peak_idx] = baseline
            continue
        
        # Use the best fit found
        amp, center, width = best_fit
        
        # Skip peaks with amplitude less than the standard deviation
        if amp <= current_std_dev:
            if verbose:
                print(f"  Fitted peak amplitude ({amp:.4f}) is below or equal to standard deviation ({std_dev:.4f}). Skipping.")
            # Mark this region to avoid repeated attempts
            residuals[peak_idx] = baseline
            continue
        
        # Store the parameters for this peak
        peak_params.extend([amp, center, width])
        
        # Calculate individual peak model
        peak_model = gaussian(t, amp, center, width)
        peak_models.append(peak_model)
        
        # Add this peak to the overall model
        overall_model += peak_model
        
        # Subtract this peak from residuals
        residuals -= peak_model
        
        # Check for convergence based on variance improvement
        current_variance = np.var(residuals)
        improvement = (prev_variance - current_variance) / initial_variance
        
        if verbose:
            print(f"  Fitted A={amp:.4f}, center={center:.6f}, width={width:.6f}, R²={best_rsq:.4f}")
            print(f"  Variance: {current_variance:.6f}, Improvement: {improvement:.2%}")
        
        # Store results without the full arrays
        results.append({
            'amplitude': amp, 
            'center': center, 
            'width': width,
            'r_squared': best_rsq,
            'improvement': improvement,
            'residual_variance': current_variance
        })
        
        # Check if residuals are approaching the noise floor
        if current_variance <= noise_floor**2:
            if verbose:
                print(f"Residuals have reached the noise floor (~{noise_floor:.3f}). Stopping.")
            break
            
        # Check if improvement is too small
        if improvement < min_improvement and iteration > 3:
            if verbose:
                print(f"Improvement below threshold ({improvement:.4%} < {min_improvement:.4%}). Stopping.")
            break
            
        prev_variance = current_variance
        
        # Additional sanity check - break if adding more peaks doesn't make sense
        if len(results) >= 3 and current_variance < 0.05 * initial_variance:
            best_fits = sorted(results, key=lambda x: x['amplitude'], reverse=True)
            smallest_amplitude = best_fits[-1]['amplitude']
            largest_amplitude = best_fits[0]['amplitude']
            # If the smallest peak is less than 10% of the largest peak, stop
            if smallest_amplitude < 0.1 * largest_amplitude and current_variance < 0.2 * initial_variance:
                if verbose:
                    print(f"Diminishing returns on peak amplitudes. Stopping.")
                break
            
        # Show progress every 10 iterations
        if verbose and (iteration + 1) % 10 == 0:
            elapsed = time.time() - start_time
            print(f"Completed {iteration+1} iterations in {elapsed:.2f} seconds")
    
    # REFINEMENT STEP: Once all peaks are identified, refine all parameters together
    if len(peak_params) > 0:
        if verbose:
            print(f"Refining all {len(peak_params)//3} peaks together for optimal fit...")
        
        try:
            # Use multi-gaussian model to refine all parameters together
            # Apply stricter bounds on centers to prevent peaks from wandering too far
            bounds_lower = []
            bounds_upper = []
            
            for i in range(len(peak_params)//3):
                idx = i * 3
                # Amplitude bounds
                bounds_lower.append(0)
                bounds_upper.append(peak_params[idx] * 2)  # Allow up to double amplitude
                
                # Center bounds - keep within reasonable window of original center
                center = peak_params[idx + 1]
                width = peak_params[idx + 2]
                bounds_lower.append(center - width)  # Allow shifting by width
                bounds_upper.append(center + width)
                
                # Width bounds
                bounds_lower.append(0.0001)
                bounds_upper.append(0.2)
            
            refined_params, _ = curve_fit(
                multi_gaussian,
                t,
                y,
                p0=peak_params,
                bounds=(bounds_lower, bounds_upper),
                maxfev=100000,
                method='trf'
            )
            
            # Update results and recalculate models
            overall_model = np.zeros_like(t)
            peak_models = []
            
            # Create temporary lists to store valid refined parameters
            valid_results = []
            valid_peak_models = []
            
            for i in range(len(refined_params) // 3):
                idx = i * 3
                amp = refined_params[idx]
                center = refined_params[idx + 1]
                width = refined_params[idx + 2]
                
                # Skip peaks with amplitude less than or equal to standard deviation
                if amp <= threshold*current_std_dev:
                    if verbose:
                        print(f"  Refined peak {i+1} amplitude ({amp:.4f}) is below or equal to standard deviation ({std_dev:.4f}). Removing from results.")
                    continue
                
                # Calculate peak model for valid peaks
                peak_model = gaussian(t, amp, center, width)
                valid_peak_models.append(peak_model)
                overall_model += peak_model
                
                # Store valid results
                peak_data = {
                    'amplitude': amp,
                    'center': center,
                    'width': width
                }
                
                # Calculate R-squared for the refined peak
                window_mask = np.abs(t - center) < window_size*2
                t_window = t[window_mask]
                y_window = y[window_mask]
                if len(t_window) >= 5:
                    y_pred = gaussian(t_window, amp, center, width)
                    ss_res = np.sum((y_window - y_pred) ** 2)
                    ss_tot = np.sum((y_window - np.mean(y_window)) ** 2)
                    peak_data['refined_r_squared'] = 1 - (ss_res / ss_tot)
                
                valid_results.append(peak_data)
            
            # Replace original lists with filtered ones
            results = valid_results
            peak_models = valid_peak_models
            
            # Recalculate final residuals
            residuals = y - overall_model
            
            if verbose:
                print(f"After refinement: {len(results)} valid peaks (amplitude > std_dev)")
                final_var = np.var(residuals)
                improvement = (initial_variance - final_var) / initial_variance
                print(f"Refinement complete. Final variance: {final_var:.6f}")
                print(f"Total improvement: {improvement:.2%}")
                
        except Exception as e:
            if verbose:
                print(f"Joint refinement failed: {e}")
                print("Using individual peak fits instead.")
    
    # Sort results by time (center)
    results.sort(key=lambda x: x['center'])
    
    peak_intervals = {
    f"Peak {i+1}": ((result['center'] - result['width'] / 2 - t[0])*24, (result['center'] + result['width'] / 2 - t[0])*24)
    for i, result in enumerate(results)
    }

    # Generate all the plots
    if verbose:
        print("Generating final plots...")
    
    # Clear the plot and redraw
    ax1.clear()
    ax2.clear()
    
    # Plot the raw data
    ax1.plot(t, y, 'k.', label='Raw Data', alpha=0.5)
    if errors is True:
        ax1.errorbar(t, y, yerr=e_y, fmt='none', ecolor='gray', capsize=4,alpha=0.4)

    ax1.axhline(threshold * current_std_dev, color='gray', ls='--', label=f'{threshold} σ')
    
    # Organize peaks by amplitude for better visualization
    peak_indices = sorted(range(len(peak_models)),key=lambda i: np.max(peak_models[i]),reverse=True)

    # Plot individual peaks (limited number for clarity)
    max_peaks_to_plot = max(8, len(peak_models))
    colors = plt.cm.viridis(np.linspace(0, 0.8, max_peaks_to_plot))
    
    if show_peaks:
        for i, idx in enumerate(peak_indices[:max_peaks_to_plot]):
                label = f'Peak {i+1}' #if i < 7 else None
                ax1.plot(t, peak_models[idx], '-', color=colors[i], alpha=0.9, label=label, linewidth=1.5, zorder=10)
                ax1.axvline(results[idx]['center'], color=colors[i], ls=':', alpha=0.7)
    
    # Plot final overall model
    ax1.plot(t, overall_model, 'r-', linewidth=2, label='Overall Model', zorder=1,alpha=0.5)
    ax1.set_ylabel('Flux [mJy]')
    ax1.legend(loc='upper right')
    ax1.set_title(f'{len(results)} Iterative Gaussian Fits')   
    ax1.grid(alpha=0.2)

    # Lower plot: Residuals
    ax2.plot(t, residuals, 'b-', label='Residuals')
    ax2.axhline(0, color='k', linestyle=':')
    ax2.set_xlabel('Time (MJD)')
    ax2.set_ylabel('Residuals')
    ax2.legend()
    ax2.grid(alpha=0.2)

    plt.tight_layout()
    plt.show()

    # Calculate overall fit quality metrics
    r_squared = 1 - (np.sum(residuals**2) / np.sum((y - np.mean(y))**2))
    rmse = np.sqrt(np.mean(residuals**2))
    
    fit_metrics = {
        'r_squared': r_squared,
        'rmse': rmse,
        'initial_variance': initial_variance,
        'final_variance': np.var(residuals),
        'improvement_percent': 100 * (1 - np.var(residuals)/initial_variance),
        'num_peaks': len(results)
    }
    
    if verbose:
        print(f"\nFit Summary:")
        print(f"Number of peaks fitted: {len(results)}")
        print(f"Overall R²: {r_squared:.4f}")
        print(f"RMSE: {rmse:.4f}")
        print(f"Variance reduction: {fit_metrics['improvement_percent']:.2f}%")
        
        print("\nPeak parameters:")
        for i, res in enumerate(results):
            print(f"Peak {i+1}: Amplitude={res['amplitude']:.2f}, "
                f"Center={res['center']:.6f}, Width={res['width']:.6f}")

    return results, fit_metrics, overall_model, peak_intervals
