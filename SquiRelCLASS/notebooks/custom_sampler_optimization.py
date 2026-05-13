import numpy as np 
import scipy 
from scipy.special import erf
import matplotlib.pyplot as plt

def lognormgaussian(q, sigma=1.0):
    """
    Returns the value of a Gaussian (normal) distribution with mean mu and standard deviation sigma at x.
    """
    coeff1 = 2.0 / np.power(2 * np.pi, 3.0)
    coeff2 = 1.0 / (q * sigma * np.sqrt(2 * np.pi))
    exponent = -0.5 * (np.log(q) / sigma) ** 2
    return coeff1 * coeff2 * np.exp(exponent)

def peak_value(sigma=1.0):
    """
    Returns the peak value of the lognormgaussian function for a given sigma.
    """
    qpeak = np.exp(-sigma**2)
    peak_value = lognormgaussian(qpeak, sigma=sigma)
    return qpeak, peak_value

def find_q_values(sigma, drop_factor=0.01):
    """
    Finds the q values where the lognormgaussian function drops to a specified percentage of its peak value.
    """
    # Find the peak value
    y_peak = peak_value(sigma)[1]

    # Calculate the target value for the drop percentage
    target_value = drop_factor * y_peak

    # Solve the quadratic equation for u = ln(q)
    a = 1 / (2 * sigma**2)
    b = 1
    c = (sigma**2 / 2) + np.log(drop_factor)

    delta = np.sqrt(b**2 - 4*a*c)
    u1 = (-b - delta) / (2*a)
    u2 = (-b + delta) / (2*a)

    q_min = np.exp(u1)
    q_max = np.exp(u2)

    return q_min, q_max, y_peak, target_value

def analytic_integral(sigma, qmin=1e-3, qmax=1e2):
    """
    Computes the analytic integral of the lognormgaussian function.
    """
    full_integral = np.exp(9 * sigma**2 / 2)/(4 * np.power(np.pi, 3))
    # The Mathematica expression needs scipy.special.erf

    # Translate the Mathematica expression to Python
    term1 = np.exp(9 * sigma**2 / 2)
    
    erf_qmax = erf((-3 * sigma**2 + np.log(qmax)) / (np.sqrt(2) * sigma))
    
    log_qmin_term = -3 * sigma**2 + np.log(qmin)
    arg_qmin = np.abs(log_qmin_term) / (np.sqrt(2) * sigma)
    erf_qmin = erf(arg_qmin) * np.sign(log_qmin_term)
    
    finite_integral = (full_integral/2) * (erf_qmax - erf_qmin)
    
    return full_integral, finite_integral

def find_drop_factor(sigma, accuracy_goal=1.0e-2, max_counter=100, drop_rescaling=0.5):
    """
    Finds the necessary drop percentage for a given accuracy goal between finite and full analytical integral.
    1. Making `accuracy_goal` smaller will force the finite integral to be closer to the full integral. This will widen the q bounds.
    2. The `drop_rescaling` parameter controls how much the drop factor is reduced each iteration. Making this smaller will slow the convergence. 
    3. The `max_counter` parameter limits the number of iterations to prevent infinite loops.
    """
    accuracy = accuracy_goal+1.0  # Initialize accuracy to a value greater than the goal
    drop_factor = 1.0e-2  # Start with a default drop factor
    counter = 0
    success = False
    while ((accuracy > accuracy_goal) and (counter < max_counter)):
        q_min, q_max, _, _ = find_q_values(sigma=sigma, drop_factor=drop_factor)
        full_integral, finite_integral = analytic_integral(sigma=sigma, qmin=q_min, qmax=q_max)
        accuracy = np.abs((finite_integral - full_integral) / full_integral)
        # print(accuracy, accuracy_goal)
        if accuracy < accuracy_goal:
            success = True
            break

        counter += 1
        drop_factor *= drop_rescaling # Making this smaller to find a more precise drop factor (i.e. lower qmax)

    if not success:
        print(f"Warning: Accuracy goal not achieved within {counter} iterations.")
        return None
    else: 
        return drop_factor

def custom_sample(f_of_q, qmin, qmax, a, n_samples=100):
    """
    Numerically samples a function f(q) between qmin and qmax using a rescaled parameter p = q^(1/a).
    This is a Monte Carlo integration.
    """
    pmin = qmin**(1/a)
    pmax = qmax**(1/a)

    # Generate uniform samples in p-space
    p_samples = np.linspace(pmin, pmax, n_samples)
    # p_samples = np.random.uniform(pmin, pmax, n_samples)  # Alternative: random uniform sampling

    # Transform samples back to q-space to evaluate f(q)
    q_samples = p_samples**a

    # The jacobian of the transformation is a * p^(a-1)
    jacobian_values = a * p_samples**(a-1)

    # Evaluate the integrand at the sample points
    integrand_values = f_of_q(q_samples) * jacobian_values

    # Values of the measure 
    dp = p_samples[1] - p_samples[0]  # Include endpoints for proper integration

    # The Monte Carlo estimate of the integral
    integral_estimate = np.sum(integrand_values * dp)

    return integral_estimate

def find_optimal_a(f_of_q, qmin, qmax, goal, tolerance, a_range=(0.1, 5.0), a_steps=50, n_max=100):
    """
    Finds the optimal 'a' that achieves the integration 'goal' within 'tolerance' 
    using the minimum number of samples.
    """
    best_a = None
    min_n_samples = float('inf')
    best_integral = None

    a_values = np.linspace(a_range[0], a_range[1], a_steps)

    for a in a_values:
        # Find the minimum n_samples for this 'a'
        n_samples = 2  # Start with a small number of samples
        while n_samples <= n_max:
            integral_estimate = custom_sample(f_of_q, qmin, qmax, a, n_samples=n_samples)
            
            if abs((integral_estimate - goal)/goal) <= tolerance:
                # This 'a' works with n_samples. Check if it's the best so far.
                if n_samples < min_n_samples:
                    min_n_samples = n_samples
                    best_a = a
                    best_integral = integral_estimate
                break  # Found min n_samples for this 'a', move to the next 'a'
            
            n_samples = int(n_samples * 1.5) # Increase n_samples for the next try

    if best_a is None:
        print(f"Warning: Could not find a suitable 'a' to meet the tolerance within n_max={n_max} samples.")

    return best_a, min_n_samples, best_integral


print("Example usage of the functions:...")
test_sigma = 1.0
test_drop_factor = find_drop_factor(sigma=test_sigma, accuracy_goal=1.0e-2)
test_qmin, test_qmax, _, _ = find_q_values(sigma=test_sigma, drop_factor=test_drop_factor)
test_full_integral, test_finite_integral = analytic_integral(sigma=test_sigma, qmin=test_qmin, qmax=test_qmax)
print("\tNecessary drop factor for accuracy goal 1.0e-2:", test_drop_factor)
print("\tq_min:", test_qmin)
print("\tq_max:", test_qmax)
print("\tFull integral:", test_full_integral)
print("\tFinite integral:", test_finite_integral)

print("Example usage of the sampling functions:...")
def f_of_q(q):
    return np.power(q, 3.0)  # Example function, can be replaced with any function of q

reference_sampling = custom_sample(f_of_q, qmin=0.01, qmax=10, a=1.0, n_samples=10000)
test_sampling = custom_sample(f_of_q, qmin=0.01, qmax=10, a=3.0, n_samples=10000)
print("\tReference sampling (a=1.0):", reference_sampling)
print("\tTest sampling (a=3.0):", test_sampling)

print("Example of finding the optimal 'a' for a given goal and tolerance...")
test_goal = test_finite_integral
test_tolerance = 1.0e-2
optimal_a, min_samples, best_integral = find_optimal_a(
    lambda q: np.power(q, 3.0) * lognormgaussian(q, sigma=test_sigma),
    qmin=test_qmin, 
    qmax=test_qmax, 
    goal=test_goal, 
    tolerance=test_tolerance
    )
print("\tOptimal 'a':", optimal_a)
print("\tMinimum samples needed:", min_samples)
print("\tBest integral estimate:", best_integral)
print("\tGoal integral:", test_goal)


# Create a plot of the lognormgaussian function
if __name__ == "__main__":
    ######################
    # Setup
    ######################
    accuracy_goal = 1e-6
    sigmas = np.arange(0.01, 2.01, 0.01)
    max_counter = 1000
    drop_rescaling = 0.9  
    a_range=(0.1, 100.0)
    a_steps=1000
    n_max=100
    results = []



    ######################
    # Build qmin, qmax lookup table. Make Plots. 
    ######################

    for sigma in sigmas:
        print(f"Processing sigma={sigma:.2f}...")
        drop_factor = find_drop_factor(
            sigma=sigma, 
            accuracy_goal=accuracy_goal, 
            max_counter=max_counter, 
            drop_rescaling=drop_rescaling
            )
        if drop_factor is not None:
            q_min, q_max, _, _ = find_q_values(sigma=sigma, drop_factor=drop_factor)
            full_integral, finite_integral = analytic_integral(sigma=sigma, qmin=q_min, qmax=q_max)
            _, finite_integral_qmin0 = analytic_integral(sigma=sigma, qmin=0., qmax=q_max)
            optimal_a, min_samples, best_integral = find_optimal_a(
                lambda q: np.power(q, 3.0) * lognormgaussian(q, sigma=sigma),
                qmin=q_min, 
                #qmin=1.e-9,  # Start from a very small qmin to avoid numerical issues
                qmax=q_max, 
                # goal=finite_integral, 
                goal=full_integral, 
                tolerance=accuracy_goal,
                a_range=a_range, 
                a_steps=a_steps,
                n_max=n_max
                )
            best_integral_qmin0 = custom_sample(
                lambda q: np.power(q, 3.0) * lognormgaussian(q, sigma=sigma),
                qmin=1.e-9,  # Consistent with find_optimal_a qmin
                qmax=q_max,
                a=optimal_a,
                n_samples=min_samples
            )
            # Only add results if find_optimal_a was successful
            if optimal_a is not None and min_samples != float('inf') and best_integral is not None:
                results.append([sigma, q_min, q_max, optimal_a, min_samples, best_integral, best_integral_qmin0, finite_integral, full_integral, np.abs((best_integral-finite_integral)/finite_integral), np.abs((best_integral-full_integral)/full_integral), np.abs((best_integral_qmin0-best_integral)/best_integral), np.abs((best_integral_qmin0-full_integral)/full_integral)])
            else:
                print(f"Warning: Could not find optimal 'a' for sigma={sigma}. Skipping.")
        else: 
            print(f"Warning: Could not find a suitable drop factor for sigma={sigma}. Skipping.")

    header = f"accuracy_goal: {accuracy_goal}\n" + \
             f"{'sigma':>22} {'qmin':>24} {'qmax':>24} {'optimal_a':>24} {'min_samples':>24} {'best_integral':>24} {'best_integral_qmin0':>24} {'finite_integral':>24} {'full_integral':>24} {'rel_err_best-finite':>24} {'rel_err_best-full':>24} {'rel_err_bestqmin0_best':>24} {'rel_err_bestqmin0_full':>24}"

    if results:
        np.savetxt("q_bounds_for_accuracy.dat", results, header=header, fmt="%.18e")
        print("Saved results to q_bounds_for_accuracy.dat")
    else:
        print("Warning: No valid results to save. Check your parameters.")

    ######################
    # Plot the lognormgaussian function for all results
    ######################    
    if results:
        # Convert results to numpy array for easier manipulation
        results_array = np.array(results)
        
        # Find the global q_min and q_max across all results
        global_q_min = np.min(results_array[:, 1])  # 2nd column is q_min
        global_q_max = np.max(results_array[:, 2])  # 3rd column is q_max
        
        plt.figure(figsize=(12, 8))
        
        # Generate q values for plotting using global bounds
        q_values = np.logspace(np.log10(global_q_min), np.log10(global_q_max), 1000)
        
        # Set up colormap
        cmap = plt.cm.magma
        n_results = len(results)
        
        # Decide which sigma values to show in legend (e.g., every 20th one, or specific indices)
        legend_indices = np.linspace(0, n_results-1, min(8, n_results), dtype=int)
        
        # Plot a line for each result
        for i, result in enumerate(results):
            sigma = result[0]
            y_values = lognormgaussian(q_values, sigma=sigma)
            color = cmap(i / (n_results - 1))  # Normalize index to [0, 1] for colormap
            
            # Only add label for selected sigma values
            label = f'σ = {sigma:.2f}' if i in legend_indices else None
            plt.loglog(q_values, y_values, color=color, label=label, linewidth=1.5)
        
        plt.xlabel('q', fontsize=14)
        plt.ylabel('lognormgaussian(q)', fontsize=14)
        plt.title('Log-Normal Gaussian Function for Different σ Values', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=10, ncol=2)  # Fewer columns since fewer legend entries
        plt.ylim(1e-20, None)
        plt.tight_layout()
        
        # Save the plot
        plt.savefig('lognormgaussian_plot.png', dpi=300, bbox_inches='tight')
        plt.savefig('lognormgaussian_plot.pdf', bbox_inches='tight')
        
        print("Plot saved as 'lognormgaussian_plot.png' and 'lognormgaussian_plot.pdf'")
        print(f"Global q range: {global_q_min:.6e} to {global_q_max:.6e}")
        
        # Optionally display the plot
        # plt.show()
        
        ######################
        # Additional plots from results
        ######################
        
        # Extract data for plotting
        sigmas = results_array[:, 0]
        q_mins = results_array[:, 1]
        q_maxs = results_array[:, 2]
        optimal_as = results_array[:, 3]
        min_samples = results_array[:, 4]
        
        # Create a 2x2 subplot figure
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Add overall title with accuracy goal
        fig.suptitle(f'Parameter Analysis (accuracy_goal = {accuracy_goal:.0e})', fontsize=16, y=0.98)
        
        # Plot 1: q_min as a function of sigma
        axes[0, 0].semilogy(sigmas, q_mins, 'b-', linewidth=2, marker='o', markersize=3)
        axes[0, 0].set_xlabel('σ', fontsize=12)
        axes[0, 0].set_ylabel('q_min', fontsize=12)
        axes[0, 0].set_title('q_min vs σ', fontsize=14)
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: q_max as a function of sigma
        axes[0, 1].semilogy(sigmas, q_maxs, 'r-', linewidth=2, marker='o', markersize=3)
        axes[0, 1].set_xlabel('σ', fontsize=12)
        axes[0, 1].set_ylabel('q_max', fontsize=12)
        axes[0, 1].set_title('q_max vs σ', fontsize=14)
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: optimal_a as a function of sigma
        axes[1, 0].plot(sigmas, optimal_as, 'g-', linewidth=2, marker='o', markersize=3)
        axes[1, 0].set_xlabel('σ', fontsize=12)
        axes[1, 0].set_ylabel('optimal_a', fontsize=12)
        axes[1, 0].set_title('Optimal a vs σ', fontsize=14)
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: min_samples as a function of sigma
        axes[1, 1].plot(sigmas, min_samples, 'm-', linewidth=2, marker='o', markersize=3)
        axes[1, 1].set_xlabel('σ', fontsize=12)
        axes[1, 1].set_ylabel('min_samples', fontsize=12)
        axes[1, 1].set_title('Min Samples vs σ', fontsize=14)
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save the parameter plots
        plt.savefig('parameter_plots.png', dpi=300, bbox_inches='tight')
        plt.savefig('parameter_plots.pdf', bbox_inches='tight')
        
        print("Parameter plots saved as 'parameter_plots.png' and 'parameter_plots.pdf'")
        
        # Optionally display the plots
        # plt.show()

        ######################
        # Fit optimal_a vs sigma with outlier removal
        ######################
        print("\nFitting optimal_a vs sigma with outlier removal...")

        # Extract data for fitting
        sigmas_fit = results_array[:, 0]
        optimal_as_fit = results_array[:, 3]

        # Outlier detection using IQR method
        Q1 = np.percentile(optimal_as_fit, 25)
        Q3 = np.percentile(optimal_as_fit, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 0.5 * IQR
        upper_bound = Q3 + 0.5 * IQR

        # Separate inliers and outliers
        is_inlier = (optimal_as_fit >= lower_bound) & (optimal_as_fit <= upper_bound)
        sigmas_inliers = sigmas_fit[is_inlier]
        optimal_as_inliers = optimal_as_fit[is_inlier]
        sigmas_outliers = sigmas_fit[~is_inlier]
        optimal_as_outliers = optimal_as_fit[~is_inlier]

        print(f"Found {len(optimal_as_outliers)} outliers.")

        # Fit a polynomial to the inliers (e.g., degree 3)
        if len(sigmas_inliers) > 3: # Need enough points to fit
            coeffs = np.polyfit(sigmas_inliers, optimal_as_inliers, 3)
            poly_fit = np.poly1d(coeffs)

            # Generate smooth x-values for plotting the fit
            sigma_fit_line = np.linspace(sigmas_fit.min(), sigmas_fit.max(), 500)
            optimal_a_fit_line = poly_fit(sigma_fit_line)

            # Create the plot
            plt.figure(figsize=(12, 8))
            
            # Plot all original data points
            plt.scatter(sigmas_fit, optimal_as_fit, color='lightblue', label='Original Data', s=30, zorder=1)
            
            # Plot the inliers
            plt.scatter(sigmas_inliers, optimal_as_inliers, color='blue', label='Inliers', s=35, zorder=2)
            
            # Plot the outliers
            if len(sigmas_outliers) > 0:
                plt.scatter(sigmas_outliers, optimal_as_outliers, color='red', marker='x', s=100, label='Outliers', zorder=3)
            
            # Plot the polynomial fit
            plt.plot(sigma_fit_line, optimal_a_fit_line, 'k--', linewidth=2, label='Polynomial Fit (deg=3)', zorder=4)
            
            # Add text showing the polynomial expression
            poly_text = f'a(σ) = {coeffs[0]:.3e}σ³ + {coeffs[1]:.3e}σ² + {coeffs[2]:.3e}σ + {coeffs[3]:.3e}'
            plt.text(0.05, 0.95, poly_text, transform=plt.gca().transAxes, fontsize=12, 
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
                     verticalalignment='top')
            
            plt.xlabel('σ', fontsize=14)
            plt.ylabel('optimal_a', fontsize=14)
            plt.title('Fit of Optimal "a" vs σ with Outlier Removal', fontsize=16)
            plt.grid(True, alpha=0.3)
            plt.legend(fontsize=12)
            plt.tight_layout()
            
            # Save the plot
            plt.savefig('optimal_a_fit_plot.png', dpi=300, bbox_inches='tight')
            plt.savefig('optimal_a_fit_plot.pdf', bbox_inches='tight')
            
            print("Fit plot saved as 'optimal_a_fit_plot.png' and 'optimal_a_fit_plot.pdf'")
            
            # Optionally display the plot
            # plt.show()

            ######################
            # Compare sampling efficiency of numerical vs. fitted 'a'
            ######################
            print("\nComparing sampling efficiency of numerical vs. fitted 'a'...")

            # Data from the original run
            sigmas_comp = results_array[:, 0]
            q_mins_comp = results_array[:, 1]
            q_maxs_comp = results_array[:, 2]
            min_samples_orig = results_array[:, 4].astype(int)
            finite_integrals_comp = results_array[:, 7]
            full_integrals_comp = results_array[:, 8]  # Use full_integral (column 8) instead of finite_integral (column 7)

            # Store the number of samples needed when using the fitted 'a'
            n_samples_fit_list = []
            extra_samples_list = []

            for i in range(len(sigmas_comp)):
                sigma = sigmas_comp[i]
                q_min = q_mins_comp[i]
                q_max = q_maxs_comp[i]
                # goal = finite_integrals_comp[i]  # Use finite_integral instead of full_integral
                goal = full_integrals_comp[i]  # Use full_integral instead of finite_integral
                n_samples_orig = min_samples_orig[i]
                
                # Get 'a' from the polynomial fit
                a_fit = poly_fit(sigma)
                
                # Now, find how many samples are needed with this a_fit to reach the tolerance
                n_samples_current = n_samples_orig # Start with the original number of samples
                
                while n_samples_current <= n_max * 5: # Allow it to search for more samples if needed
                    integral_estimate = custom_sample(
                        lambda q: np.power(q, 3.0) * lognormgaussian(q, sigma=sigma),
                        qmin=q_min,
                        qmax=q_max,
                        a=a_fit,
                        n_samples=n_samples_current
                    )
                    
                    if abs((integral_estimate - goal) / goal) <= accuracy_goal:
                        # Tolerance met
                        break
                    
                    # Increase samples if tolerance not met
                    n_samples_current = int(n_samples_current * 1.5) if n_samples_current > 0 else 2
                
                if n_samples_current > n_max * 5:
                    print(f"Warning: Could not meet tolerance with fitted 'a' for sigma={sigma:.2f}. Using NaN.")
                    n_samples_fit_list.append(np.nan)
                    extra_samples_list.append(np.nan)
                else:
                    n_samples_fit_list.append(n_samples_current)
                    extra_samples_list.append(n_samples_current - n_samples_orig)

            # Plot the comparison of sampling efficiency
            fig, axes = plt.subplots(2, 1, figsize=(12, 12), sharex=True)
            fig.suptitle('Comparison of Sampling Efficiency', fontsize=18, y=0.95)

            # Plot 1: Extra samples needed
            axes[0].plot(sigmas_comp, extra_samples_list, 'o-', label='Extra Samples (Fit - Numerical)')
            axes[0].axhline(0, color='k', linestyle='--', linewidth=1) # Add a line at y=0 for reference
            axes[0].set_ylabel('Extra Samples (Fit - Numerical)', fontsize=14)
            axes[0].set_title('Relative Cost', fontsize=16)
            axes[0].grid(True, alpha=0.3)
            axes[0].legend(fontsize=12)

            # Plot 2: Absolute number of samples
            axes[1].plot(sigmas_comp, min_samples_orig, 'o-', color='g', label='Samples with Numerical "a"')
            axes[1].plot(sigmas_comp, n_samples_fit_list, 'o-', color='purple', label='Samples with Fitted "a"')
            axes[1].set_xlabel('σ', fontsize=14)
            axes[1].set_ylabel('Total Sample Points', fontsize=14)
            axes[1].set_title('Absolute Cost', fontsize=16)
            axes[1].set_yscale('log')
            axes[1].grid(True, alpha=0.3)
            axes[1].legend(fontsize=12)
            
            plt.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust layout to make room for suptitle

            # Save the plot
            plt.savefig('sampling_efficiency_comparison.png', dpi=300, bbox_inches='tight')
            plt.savefig('sampling_efficiency_comparison.pdf', bbox_inches='tight')

            print("Sampling efficiency comparison plot saved as 'sampling_efficiency_comparison.png' and 'sampling_efficiency_comparison.pdf'")

            # Optionally display the plot
            # plt.show()
        else:
            print("Not enough inliers to perform a fit.")
        
    else:
        print("No results to plot.")

    # ######################
    # # Generate 2D heatmap: min_samples vs sigma and accuracy_goal
    # ######################
    
    # print("\nGenerating 2D heatmap data...")
    
    # # Define parameter grids
    # sigma_grid = np.linspace(0.01, 2.0, 40)  # 40 points for sigma
    # accuracy_grid = np.logspace(np.log10(0.1), np.log10(0.00000001), 24)  # 20 points for accuracy_goal
    
    # # Initialize the result grid
    # min_samples_grid = np.full((len(accuracy_grid), len(sigma_grid)), np.nan)
    
    # # Calculate min_samples for each combination
    # total_combinations = len(sigma_grid) * len(accuracy_grid)
    # current_combination = 0
    
    # for i, accuracy in enumerate(accuracy_grid):
    #     for j, sigma in enumerate(sigma_grid):
    #         current_combination += 1
    #         print(f"Processing combination {current_combination}/{total_combinations}: σ={sigma:.3f}, accuracy={accuracy:.2e}")
            
    #         # Find drop factor for this accuracy goal
    #         drop_factor = find_drop_factor(
    #             sigma=sigma, 
    #             accuracy_goal=accuracy, 
    #             max_counter=max_counter, 
    #             drop_rescaling=drop_rescaling
    #         )
            
    #         if drop_factor is not None:
    #             q_min, q_max, _, _ = find_q_values(sigma=sigma, drop_factor=drop_factor)
    #             full_integral, finite_integral = analytic_integral(sigma=sigma, qmin=q_min, qmax=q_max)
                
    #             optimal_a, min_samples, best_integral = find_optimal_a(
    #                 lambda q: np.power(q, 3.0) * lognormgaussian(q, sigma=sigma),
    #                 qmin=q_min, 
    #                 qmax=q_max, 
    #                 goal=finite_integral, 
    #                 tolerance=accuracy,
    #                 a_range=a_range, 
    #                 a_steps=a_steps,
    #                 n_max=n_max
    #             )
                
    #             if optimal_a is not None and min_samples != float('inf'):
    #                 min_samples_grid[i, j] = min_samples
    #             else:
    #                 print(f"    Warning: Could not find optimal 'a' for σ={sigma:.3f}, accuracy={accuracy:.2e}")
    #         else:
    #             print(f"    Warning: Could not find drop factor for σ={sigma:.3f}, accuracy={accuracy:.2e}")
    
    # # Save the heatmap data
    # heatmap_data = {
    #     'sigma_grid': sigma_grid,
    #     'accuracy_grid': accuracy_grid,
    #     'min_samples_grid': min_samples_grid
    # }
    # np.savez('heatmap_data.npz', **heatmap_data)
    # print("Heatmap data saved as 'heatmap_data.npz'")
    
    # # Create the heatmap plot
    # plt.figure(figsize=(12, 8))
    
    # # Create meshgrid for plotting
    # Sigma, Accuracy = np.meshgrid(sigma_grid, accuracy_grid)
    
    # # Create the heatmap
    # im = plt.pcolormesh(Sigma, Accuracy, min_samples_grid, cmap='viridis', shading='auto')
    
    # # Add colorbar
    # cbar = plt.colorbar(im)
    # cbar.set_label('min_samples', fontsize=14)
    
    # # Set labels and title
    # plt.xlabel('σ', fontsize=14)
    # plt.ylabel('accuracy_goal', fontsize=14)
    # plt.title('Min Samples Required vs σ and Accuracy Goal', fontsize=16)
    
    # # Use log scale for y-axis since accuracy_goal spans many orders of magnitude
    # plt.yscale('log')
    
    # # Add grid
    # plt.grid(True, alpha=0.3)
    
    # plt.tight_layout()
    
    # # Save the heatmap
    # plt.savefig('min_samples_heatmap.png', dpi=300, bbox_inches='tight')
    # plt.savefig('min_samples_heatmap.pdf', bbox_inches='tight')
    
    # print("Heatmap saved as 'min_samples_heatmap.png' and 'min_samples_heatmap.pdf'")
    
    # # Optionally display the plot
    # # plt.show()

