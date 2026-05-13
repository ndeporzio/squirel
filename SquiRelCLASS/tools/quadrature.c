/******************************************/
/* Quadrature Sampling Strategy for CLASS */
/* 10/12 2010                             */
/* Thomas Tram                            */
/******************************************/
#include "quadrature.h"
#include "background.h"

int get_qsampling_manual(double *x,
			 double *w,
			 int N,
			 double qmax,
       double a, 
			 enum ncdm_quadrature_method method,
			 double *qvec,
			 int qsiz,
			 int (*function)(void * params_for_function, double q, double *f0),
			 void * params_for_function,
			 ErrorMsg errmsg) {

  double y, h, t;
  double *b, *c;
  int i;

  switch (method){ 
  case (qm_auto) :
    return _FAILURE_;
  case (qm_Laguerre) :
    /* Allocate storage for Laguerre coefficients: */
    class_alloc(b,N*sizeof(double),errmsg);
    class_alloc(c,N*sizeof(double),errmsg);
    compute_Laguerre(x,w,N,0.0,b,c,_TRUE_);
    for (i=0; i<N; i++){
      (*function)(params_for_function,x[i],&y);
      w[i] *= y;
    }
    free(b);
    free(c);
    return _SUCCESS_;
  case (qm_trapz) :
    for (i=0; i<N; i++){
      /** Note that we count q=0 as an extra point with weight 0 */
      h = qmax/N;
      x[i] = h + i*h;
      (*function)(params_for_function,x[i],&y);
      w[i] = y*h;
      if (i==N-1)
	w[i] *=0.5;
    }
    return _SUCCESS_;
  case (qm_trapz_indefinite) :
    /** We do the variable transformation q = 1/t-1. The trapezoidal rule is closed, but since the distribution function
	goes to zero in both limits, we can use an effectively N+2 rule simply by not using the exterior points. */
    for (i=0; i<N; i++){
      h = 1.0/(N+1.0);
      t = h + i*h;
      x[i] = 1.0/t-1.0;
      (*function)(params_for_function,x[i],&y);
      w[i] = y*h/t/t;
    }
    return _SUCCESS_;
  case (qm_trapz_indefinite_scaled) :
    /** We do the variable transformation q = 1/t-1. The trapezoidal rule is closed, but since the distribution function
	goes to zero in both limits, we can use an effectively N+2 rule simply by not using the exterior points. */
    
    /* Check for log-normal distribution optimization before the loop */
    double optimal_qmax = qmax;
    double optimal_qmin = 0.0;  /* Will be set from lookup table for log-normal distributions */
    int optimal_N = N;  /* Start with user value, will be optimized if log-normal detected */
    double optimal_a = a;
    
    struct background_parameters_for_distributions * pbadist_local = 
        (struct background_parameters_for_distributions *) params_for_function;
    struct background * pba_local = pbadist_local->pba;
    double *param = pba_local->ncdm_psd_parameters;
    
    /* If we have parameters and this is a log-normal distribution (param[0] == 2) */
    if (param != NULL && param[0] == 2.0 && pbadist_local->n_ncdm == 0) {
      double sigma = param[1];
      
      /* Polynomial fit coefficients for optimal 'a' from custom_sampler_optimization.py */
      double c0 = -1.397e+01;  /* σ³ coefficient */
      double c1 = 4.018e+01;   /* σ² coefficient */  
      double c2 = -3.187e+01;  /* σ coefficient */
      double c3 = 1.079e+01;   /* constant term */
      
      /* Calculate optimal 'a' using polynomial fit: a(σ) = c₀σ³ + c₁σ² + c₂σ + c₃ */
      double polynomial_a = c0*sigma*sigma*sigma + c1*sigma*sigma + c2*sigma + c3;
      optimal_a = polynomial_a;
      
      /* Ensure reasonable bounds for 'a' */
      if (optimal_a < 0.01) optimal_a = 0.01;
      if (optimal_a > 20.0) optimal_a = 20.0;  /* Increased upper bound for larger sigma values */

      /* Lookup table from q_bounds_for_accuracy.dat (sigma, qmin, qmax, optimal_a, min_samples) */
      /* Data covers sigma from 0.01 to 0.17 with accuracy goal 1e-6 */
      /* Using 20 sample points for all sigma values - ensures at least 1e-6 tolerance */
      double sigma_table[] = {
        0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 , 0.11,
        0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2 , 0.21, 0.22,
        0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3 , 0.31, 0.32, 0.33,
        0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4 , 0.41, 0.42, 0.43, 0.44,
        0.45, 0.46, 0.47, 0.48, 0.49, 0.5 , 0.51, 0.52, 0.53, 0.54, 0.55,
        0.56, 0.57, 0.58, 0.59, 0.6 , 0.61, 0.62, 0.63, 0.64, 0.65, 0.66,
        0.67, 0.68, 0.69, 0.7 , 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77,
        0.78, 0.79, 0.8 , 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88,
        0.89, 0.9 , 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
        1.  , 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1 ,
        1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.2 , 1.21,
        1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3 , 1.31, 1.32,
        1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4 , 1.41, 1.42, 1.43,
        1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5 , 1.51, 1.52, 1.53, 1.54,
        1.55, 1.56, 1.57, 1.58, 1.59, 1.6 , 1.61, 1.62, 1.63, 1.64, 1.65,
        1.66, 1.67, 1.68, 1.69, 1.7 , 1.71, 1.72, 1.73, 1.74, 1.75, 1.76,
        1.77, 1.78, 1.79, 1.8 , 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87,
        1.88, 1.89, 1.9 , 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98,
        1.99, 2.
      };
      double qmin_table[] = {
        9.51929260e-01, 9.05988100e-01, 8.61538812e-01, 8.18759228e-01,
        7.77620589e-01, 7.37161850e-01, 6.98097031e-01, 6.61529564e-01,
        6.25336624e-01, 5.90537449e-01, 5.57126630e-01, 5.25094042e-01,
        4.94425263e-01, 4.65101979e-01, 4.37102380e-01, 4.10401543e-01,
        3.84971795e-01, 3.60783068e-01, 3.37803232e-01, 3.15998408e-01,
        2.95333271e-01, 2.75771331e-01, 2.56180269e-01, 2.38748808e-01,
        2.22307941e-01, 2.06818576e-01, 1.92241624e-01, 1.77645013e-01,
        1.64816488e-01, 1.52784564e-01, 1.41511508e-01, 1.30231501e-01,
        1.20403862e-01, 1.11225758e-01, 1.02662772e-01, 9.41037286e-02,
        8.67055290e-02, 7.93172223e-02, 7.29541227e-02, 6.70479636e-02,
        6.11566825e-02, 5.61094073e-02, 5.10798716e-02, 4.67849866e-02,
        4.25096268e-02, 3.88702649e-02, 3.52513163e-02, 3.21800877e-02,
        2.91294344e-02, 2.65481005e-02, 2.39869366e-02, 2.18259533e-02,
        1.96843093e-02, 1.77358446e-02, 1.60984553e-02, 1.44788786e-02,
        1.30100056e-02, 1.17804041e-02, 1.05665781e-02, 9.46906415e-03,
        8.55373224e-03, 7.65194488e-03, 6.83900335e-03, 6.10689889e-03,
        5.49917617e-03, 4.90212126e-03, 4.36599485e-03, 3.88505445e-03,
        3.45404232e-03, 3.06815115e-03, 2.74992074e-03, 2.43864511e-03,
        2.16073144e-03, 1.91283940e-03, 1.69193428e-03, 1.49526250e-03,
        1.32032861e-03, 1.16487390e-03, 1.02685637e-03, 9.04432137e-04,
        7.95938167e-04, 6.99876279e-04, 6.14898352e-04, 5.39792672e-04,
        4.73471350e-04, 4.14958755e-04, 3.63380889e-04, 3.14414361e-04,
        2.74867129e-04, 2.40100734e-04, 2.09563550e-04, 1.82764103e-04,
        1.57435130e-04, 1.37073855e-04, 1.19251677e-04, 1.03665080e-04,
        8.89861339e-05, 7.72287651e-05, 6.69726887e-05, 5.80335895e-05,
        4.96445455e-05, 4.29489253e-05, 3.71278392e-05, 3.16794112e-05,
        2.73420431e-05, 2.32897232e-05, 2.00691766e-05, 1.72809275e-05,
        1.46826861e-05, 1.26229627e-05, 1.07070204e-05, 9.19067180e-06,
        7.78266662e-06, 6.67013058e-06, 5.63889733e-06, 4.76304034e-06,
        4.07262988e-06, 3.43439960e-06, 2.93210489e-06, 2.46857160e-06,
        2.07657743e-06, 1.76879477e-06, 1.48551444e-06, 1.24656777e-06,
        1.05938780e-06, 8.87559467e-07, 7.42989751e-07, 6.29999848e-07,
        5.26545003e-07, 4.39720374e-07, 3.66914418e-07, 3.05915234e-07,
        2.58418676e-07, 2.15119561e-07, 1.78931597e-07, 1.48712098e-07,
        1.23497645e-07, 1.02476700e-07, 8.49662871e-08, 7.14016780e-08,
        5.91102399e-08, 4.88961698e-08, 4.04152896e-08, 3.33792212e-08,
        2.75465525e-08, 2.27153669e-08, 1.87169286e-08, 1.54103541e-08,
        1.26781159e-08, 1.02702638e-08, 8.43597193e-09, 6.92397232e-09,
        5.67862277e-09, 4.65370701e-09, 3.81087050e-09, 3.11830998e-09,
        2.51175252e-09, 2.05208476e-09, 1.67527500e-09, 1.36662720e-09,
        1.11400700e-09, 8.93724436e-10, 7.27399339e-10, 5.91586682e-10,
        4.73471206e-10, 3.84481241e-10, 3.11986067e-10, 2.49101901e-10,
        2.01826291e-10, 1.63402645e-10, 1.30159112e-10, 1.05220657e-10,
        8.49982943e-11, 6.75470163e-11, 5.44834722e-11, 4.32295215e-11,
        3.48169948e-11, 2.75822409e-11, 2.21817314e-11, 1.75452352e-11,
        1.40890877e-11, 1.11269193e-11, 8.78063218e-12, 7.03529420e-12,
        5.54326441e-12, 4.43492780e-12, 3.48903307e-12, 2.74274750e-12,
        2.18952302e-12, 1.71857887e-12, 1.34788745e-12, 1.07365767e-12,
        8.40800799e-13, 6.57940416e-13, 5.14454325e-13, 4.08595261e-13,
        3.19008165e-13, 2.48873685e-13, 1.94010612e-13, 1.51126959e-13
      }; 
      double qmax_table[] = {
        1.05028815e+00, 1.10288460e+00, 1.15862641e+00, 1.21745817e+00,
        1.27956036e+00, 1.34682208e+00, 1.41849603e+00, 1.49242245e+00,
        1.57344137e+00, 1.65984168e+00, 1.75200829e+00, 1.85035573e+00,
        1.95533058e+00, 2.06741408e+00, 2.18712486e+00, 2.31502208e+00,
        2.45170869e+00, 2.59783504e+00, 2.75410275e+00, 2.92126898e+00,
        3.10015102e+00, 3.29163127e+00, 3.51160749e+00, 3.73274278e+00,
        3.96970481e+00, 4.22370757e+00, 4.49606351e+00, 4.81226577e+00,
        5.12803537e+00, 5.46698037e+00, 5.83091675e+00, 6.25662960e+00,
        6.67990437e+00, 7.13486467e+00, 7.62403471e+00, 8.20019234e+00,
        8.77088423e+00, 9.44513690e+00, 1.01120167e+01, 1.08302922e+01,
        1.16827821e+01, 1.25240625e+01, 1.35253373e+01, 1.45122472e+01,
        1.56900180e+01, 1.68495798e+01, 1.82369984e+01, 1.96015258e+01,
        2.12382993e+01, 2.28464805e+01, 2.47802238e+01, 2.66784333e+01,
        2.89663104e+01, 3.14678872e+01, 3.39209207e+01, 3.68871901e+01,
        4.01345097e+01, 4.33158146e+01, 4.71747726e+01, 5.14044734e+01,
        5.55446706e+01, 6.05818553e+01, 6.61095534e+01, 7.21780643e+01,
        7.81130383e+01, 8.53608237e+01, 9.33269763e+01, 1.02086274e+02,
        1.11721598e+02, 1.22324840e+02, 1.32685702e+02, 1.45403506e+02,
        1.59414620e+02, 1.74856658e+02, 1.91882433e+02, 2.10661695e+02,
        2.31383065e+02, 2.54256205e+02, 2.79514232e+02, 3.07416432e+02,
        3.38251281e+02, 3.72339839e+02, 4.10039536e+02, 4.51748420e+02,
        4.97909908e+02, 5.49018111e+02, 6.05623794e+02, 6.75868697e+02,
        7.46221702e+02, 8.24231963e+02, 9.10763740e+02, 1.00678191e+03,
        1.12630434e+03, 1.24612086e+03, 1.37922133e+03, 1.52712970e+03,
        1.71167893e+03, 1.89682851e+03, 2.10280123e+03, 2.33201641e+03,
        2.61867052e+03, 2.90648269e+03, 3.22711390e+03, 3.62875023e+03,
        4.03227092e+03, 4.53825728e+03, 5.04687956e+03, 5.61449609e+03,
        6.32738032e+03, 7.04443322e+03, 7.94589660e+03, 8.85308729e+03,
        9.99471293e+03, 1.11441628e+04, 1.25920636e+04, 1.42345070e+04,
        1.58896241e+04, 1.79773594e+04, 2.00822175e+04, 2.27397750e+04,
        2.57602913e+04, 2.88081083e+04, 3.26614062e+04, 3.70460335e+04,
        4.14738904e+04, 4.70793762e+04, 5.34650574e+04, 5.99188400e+04,
        6.80997644e+04, 7.74297867e+04, 8.80744017e+04, 1.00223481e+05,
        1.12519834e+05, 1.28139125e+05, 1.45985450e+05, 1.66383996e+05,
        1.89708460e+05, 2.16388434e+05, 2.46917927e+05, 2.77879950e+05,
        3.17318621e+05, 3.62495294e+05, 4.14263560e+05, 4.73606578e+05,
        5.41657111e+05, 6.19720710e+05, 7.09302554e+05, 8.12138533e+05,
        9.30231245e+05, 1.08166614e+06, 1.23992432e+06, 1.42185967e+06,
        1.63108681e+06, 1.87178247e+06, 2.14877441e+06, 2.46764461e+06,
        2.87765162e+06, 3.30719053e+06, 3.80219278e+06, 4.37282596e+06,
        5.03086458e+06, 5.87856725e+06, 6.76819839e+06, 7.79515463e+06,
        9.11954086e+06, 1.05109288e+07, 1.21187280e+07, 1.41944140e+07,
        1.63774137e+07, 1.89024811e+07, 2.21658522e+07, 2.56014993e+07,
        2.95794354e+07, 3.47259298e+07, 4.01496276e+07, 4.71713600e+07,
        5.45765844e+07, 6.41701685e+07, 7.42948414e+07, 8.74203544e+07,
        1.01282130e+08, 1.19264441e+08, 1.40493630e+08, 1.62938210e+08,
        1.92082888e+08, 2.22916937e+08, 2.62982087e+08, 3.10365428e+08,
        3.60548627e+08, 4.25818350e+08, 5.03091491e+08, 5.85016093e+08,
        6.91671455e+08, 8.18072908e+08, 9.67929313e+08, 1.12702957e+09,
        1.33441616e+09, 1.58053766e+09, 1.87273082e+09, 2.21974047e+09
      };
      int N_table[] = {
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 
        20, 20, 20, 20, 20, 20, 20, 20
      };
      int table_size = 200;
      
      /* Linear interpolation for optimal values */
      if (sigma >= sigma_table[0] && sigma <= sigma_table[table_size-1]) {
        /* Find bracketing indices */
        int idx = 0;
        while (idx < table_size-1 && sigma > sigma_table[idx+1]) {
          idx++;
        }
        
        if (idx < table_size-1) {
          /* Linear interpolation */
          double t = (sigma - sigma_table[idx]) / (sigma_table[idx+1] - sigma_table[idx]);
          optimal_qmax = qmax_table[idx] + t * (qmax_table[idx+1] - qmax_table[idx]);
          optimal_qmin = qmin_table[idx] + t * (qmin_table[idx+1] - qmin_table[idx]);
          optimal_N = (t < 0.5) ? N_table[idx] : N_table[idx+1];  /* Use nearest neighbor for N */
        } else {
          /* Use last values */
          optimal_qmax = qmax_table[table_size-1];
          optimal_qmin = qmin_table[table_size-1];
          optimal_N = N_table[table_size-1];
        }
      } else if (sigma < sigma_table[0]) {
        /* Use first values for sigma below range */
        optimal_qmax = qmax_table[0];
        optimal_qmin = qmin_table[0];
        optimal_N = N_table[0];
      } else {
        /* Use last values for sigma above range */
        optimal_qmax = qmax_table[table_size-1];
        optimal_qmin = qmin_table[table_size-1];
        optimal_N = N_table[table_size-1];
      }
      
      /* Use the input N value directly - let user control via ncdm_N_momentum_bins */
      optimal_N = N;
      
      /* Print info about optimization */
      printf("Log-normal distribution detected: σ=%.3f\n", sigma);
      printf("  Polynomial-optimized 'a': %.3f (polynomial predicted: %.3f, input 'a': %.3f)\n", 
             optimal_a, polynomial_a, a);
      printf("  Lookup-optimized qmin: %.6e (lookup table)\n", optimal_qmin);
      printf("  Lookup-optimized qmax: %.3f (input qmax: %.3f)\n", 
             optimal_qmax, qmax);
      printf("  Optimized q-range: [%.6e, %.3f] vs input: [0.0, %.3f]\n", 
             optimal_qmin, optimal_qmax, qmax);
      printf("  Using N=%d momentum bins (set via ncdm_N_momentum_bins parameter)\n", optimal_N);
    } else {
      /* For non-log-normal distributions, use traditional sampling from 0 */
      optimal_qmin = 0.0;
    }
    
    for (i=0; i<optimal_N; i++){
      /** ND 5-20-25: This method rescales the sampling in q by an arbitrary polynomial sampling function. 
       * (e.g. q^2, q^3, etc.). This is useful for sampling the distribution function more densely in the low q region.
       * Particularly, by careful choice of qmax, the number of samples, and the scaling function, you can sample 
       * distributions densely at low q (where the distribution is changing significantly) and sparsely at high q (where 
       * the distribution is not changing significantly, but still has a significant amplitude), for example. 
       * 
       * ENHANCED 7-23-25: For log-normal distributions, automatically use optimized 'a' parameter based on sigma.
       * ENHANCED 8-5-25: For log-normal distributions, use optimized qmin and qmax from lookup table.
       */
       
      /* For log-normal optimization, sample between optimal_qmin and optimal_qmax */
      if (param != NULL && param[0] == 2.0 && pbadist_local->n_ncdm == 0) {
        /* Transform to sample between qmin and qmax using the scaling parameter 'a' */
        double t_min = pow(optimal_qmin, 1.0/optimal_a);
        double t_max = pow(optimal_qmax, 1.0/optimal_a);
        double dt = (t_max - t_min) / optimal_N;
        t = t_min + (i+1) * dt;
        x[i] = pow(t, optimal_a);
        (*function)(params_for_function,x[i],&y);
        /* Jacobian includes the transformation and the sampling density */
        w[i] = y * optimal_a * pow(t, optimal_a-1) * dt;
      } else {
        /* Original implementation for non-log-normal distributions */
        double dt = pow(optimal_qmax, 1.0/optimal_a) / optimal_N;
        t = (i+1) * dt;
        x[i] = pow(t, optimal_a);
        (*function)(params_for_function,x[i],&y);
        w[i] = y * optimal_a * pow(t, optimal_a-1) * dt;
      }
      /** End modification */
      /**printf("(t,q,y,w) = (%g,%g%g,%g)\n",t,x[i],y,w[i]);*/
    }
    
    /* Fill remaining array elements with zeros if we used fewer points than allocated */
    for (i=optimal_N; i<N; i++){
      x[i] = 0.0;
      w[i] = 0.0;
    }
    return _SUCCESS_;
  case (qm_lognormal_optimized) :
    /** Optimized sampling for log-normal distributions using fitted polynomial for 'a' parameter.
     * This method implements the optimal 'a' value based on sigma parameter from log-normal distribution.
     * The polynomial fit: a(σ) = c₀σ³ + c₁σ² + c₂σ + c₃ was derived from numerical optimization.
     */
    for (i=0; i<N; i++){
      /* Extract sigma parameter from ncdm_psd_parameters if available */
      struct background_parameters_for_distributions * pbadist_local = 
          (struct background_parameters_for_distributions *) params_for_function;
      struct background * pba_local = pbadist_local->pba;
      double *param = pba_local->ncdm_psd_parameters;
      double sigma = 1.0; /* default value */
      double optimal_a = a; /* default to user-provided value */
      
      /* If we have parameters and this is a log-normal distribution (param[0] == 2) */
      if (param != NULL && param[0] == 2.0 && pbadist_local->n_ncdm == 0) {
        sigma = param[1];
        /* Polynomial fit coefficients from optimization */
        double c0 = -1.397e+01;  /* σ³ coefficient */
        double c1 = 4.018e+01;   /* σ² coefficient */  
        double c2 = -3.187e+01;  /* σ coefficient */
        double c3 = 1.079e+01;   /* constant term */
        
        /* Calculate optimal 'a' using polynomial fit */
        optimal_a = c0*sigma*sigma*sigma + c1*sigma*sigma + c2*sigma + c3;
        
        /* Ensure reasonable bounds for 'a' */
        if (optimal_a < 0.1) optimal_a = 0.1;
        if (optimal_a > 100.0) optimal_a = 100.0;
      }
      
      /* Use the optimal (or default) 'a' value for sampling */
      double dt = pow(qmax, 1.0/optimal_a) / N;
      t = (i+1) * dt;
      x[i] = pow(t, optimal_a);
      (*function)(params_for_function,x[i],&y);
      w[i] = y * optimal_a * pow(t, optimal_a-1) * dt;
    }
    return _SUCCESS_;
  }
  return _SUCCESS_;
}
int get_qsampling(double *x,
		  double *w,
		  int *N,
		  int N_max,
		  double rtol,
		  double *qvec,
		  int qsiz,
		  int (*test)(void * params_for_function, double q, double *psi),
		  int (*function)(void * params_for_function, double q, double *f0),
		  void * params_for_function,
		  ErrorMsg errmsg) {

  /* This routine returns the fewest possible number of abscissas and weights under
     the requirement that a test function folded with the neutrino distribution function
     can be integrated to an accuracy of rtol. If the distribution function is Fermi-Dirac
     or close, a Laguerre quadrature formula is often the best choice.

     This function combines two completely different strategies: Adaptive Gauss-Kronrod
     quadrature and Laguerres quadrature formula. */

  int i, NL=2,NR,level,Nadapt=0,NLag,NLag_max,Nold=NL;
  int adapt_converging=_FALSE_,Laguerre_converging=_FALSE_,combined_converging=_FALSE_;
  double y,y1,y2,I,Igk,err,ILag,*b,*c;
  qss_node *root,*root_comb;
  double I_comb,I_atzero,I_atinf,I_comb2;
  int N_comb=0,N_comb_lag=16,N_comb_leg=4;
  double a_comb,b_comb;
  double q_leg[4],w_leg[4];
  double q_lag[N_comb_lag],w_lag[N_comb_lag];
  char method_chosen[40];
  double qmin=0., qmax=0., qmaxm1=0.;
  double *wcomb2=NULL,delq;
  double Itot=0.0;
  int zeroskip=0;

  /* Set roots and weights for Gauss-Legendre 4 point rule: */
  q_leg[3] = sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0);
  q_leg[2] = sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0);
  q_leg[1] = -q_leg[2];
  q_leg[0] = -q_leg[3];
  w_leg[3] = (18.0-sqrt(30.0))/36.0;
  w_leg[2] = (18.0+sqrt(30.0))/36.0;
  w_leg[1] = w_leg[2];
  w_leg[0] = w_leg[3];

  /* Allocate storage for Laguerre coefficients: */
  class_alloc(b,N_max*sizeof(double),errmsg);
  class_alloc(c,N_max*sizeof(double),errmsg);
  /* If a vector of q values has been passed, use it: */
  if ((qvec!=NULL)&&(qsiz>1)){
    qmin = qvec[0];
    qmax = qvec[qsiz-1];
    qmaxm1 = qvec[qsiz-2];
  }
  else{
    qvec = NULL;
  }

  /* First do the adaptive quadrature - this will also give the value of the integral: */
  gk_adapt(&root,(*test),(*function), params_for_function,
	   rtol*1e-4, 1, 0.0, 1.0, _TRUE_, errmsg);
  /* Do a leaf count: */
  leaf_count(root);
  /* I can get the integral now: */
  I = get_integral(root, 1);
  //printf("I = %le |%le, used points: %d\n",I,I-1.0,15*root->leaf_childs);
  Itot += I;

  /* Starting from the top, move down in levels until tolerance is met: */
  for(level=root->leaf_childs; level>=1; level--){
    Igk = get_integral(root,level);
    err = I-Igk;
    if (fabs(err/Itot)<rtol) break;
  }
  if (level>0){
    /* Reduce tree to the found level:*/
    reduce_tree(root,level);
    /* Count the new leafs: */
    leaf_count(root);
    /* I know know how many function evaluations is
       required by the adaptively sampled grid:*/
    Nadapt = 15*root->leaf_childs;
    /* The adaptive routine could not recieve required precision
       using less than the required maximal number of points.*/
    if (Nadapt <= N_max) adapt_converging = _TRUE_;
  }


  /* Combined adaptive quadrature and Laguerre rescaled quadrature: */
  if(qvec!=NULL){
    /* Evaluate [0;qmin] using 4 point Gauss-Legendre rule: */
    (*function)(params_for_function,qmin,&y2);
    for(i=0,I_atzero=0.0; i<N_comb_leg; i++){
      q_leg[i] = 0.5*qmin*(1.0+q_leg[i]);
      w_leg[i] = w_leg[i]*0.5*qmin*y2;
      (*test)(params_for_function,q_leg[i],&y);
      I_atzero +=w_leg[i]*y;
    }

    /* Find asymptotic extrapolation:*/
    (*function)(params_for_function,qmaxm1,&y1);
    (*function)(params_for_function,qmax,&y2);

    b_comb = (y1/y2-1.0)/(qmax-qmaxm1);
    b_comb = MAX(b_comb,1e-100);
    //c_comb = -b_comb*qmax;
    a_comb = y2*exp(b_comb*qmax);
    // printf("f(q) = %g*exp(-%g*q) \n",a_comb,b_comb);
    //(*function)(params_for_function,100,&y2);
    //printf("f(100) = %e ?= %e\n",y2,a_comb*exp(-b_comb*100));

    /* Evaluate tail using 6 point Laguerre: */
    compute_Laguerre(q_lag,w_lag,N_comb_lag,0.0,b,c,_TRUE_);
    for (i=0,I_atinf=0.0; i<N_comb_lag; i++){
      w_lag[i] *= exp(-q_lag[i]);
      q_lag[i] = qmax + q_lag[i]/b_comb;
      w_lag[i] = a_comb/b_comb*exp(-b_comb*qmax)*w_lag[i];
      (*test)(params_for_function,q_lag[i],&y);
      I_atinf +=w_lag[i]*y;
    }

    /* Do the adaptive quadrature - this will also give the main part of the integral: */
    gk_adapt(&root_comb,(*test),(*function), params_for_function,
	     rtol*1e-2, 1, qmin, qmax, _FALSE_, errmsg);
    /* Do a leaf count: */
    leaf_count(root_comb);
    /* Starting from the top, move down in levels until tolerance is met: */
    for(level=root_comb->leaf_childs; level>=1; level--){
      I_comb = get_integral(root_comb,level);
      //printf("%le + %le + %le = %le | %le\n",
      //     I_atzero,I_atinf,I_comb,I_comb+I_atinf+I_atzero,I_comb+I_atinf+I_atzero-1.0);
      I_comb +=(I_atinf+I_atzero);
      err = I-I_comb;
      if (fabs(err/Itot)<rtol) break;
    }
    /* Reduce tree to the found level:*/
    if (level>0){
      reduce_tree(root_comb,level);
      /* Count the new leafs: */
      leaf_count(root_comb);
      N_comb = 15*root_comb->leaf_childs+N_comb_lag+N_comb_leg;
      if (N_comb <= N_max) combined_converging = _TRUE_;
    }

    /* Do the second combined quadrature: Same as above, but with trapezoidal rule
       using the given q grid:  */
    class_alloc(wcomb2,qsiz*sizeof(double),errmsg);
    I_comb2 = 0.0;
    for(i=0; i<qsiz; i++){
      if(i==0){
	delq = qvec[1]-qvec[0];
      }
      else if(i==qsiz-1){
	delq = qvec[qsiz-1]-qvec[qsiz-2];
      }
      else{
	delq = qvec[i+1]-qvec[i-1];
      }
      (*function)(params_for_function,qvec[i],&y2);
      wcomb2[i] = 0.5*y2*delq;
      (*test)(params_for_function,qvec[i],&y);
      I_comb2 +=wcomb2[i]*y;
    }
    I_comb2 +=(I_atzero+I_atinf);
    err = I - I_comb2;
    //    if(fabs(err/Itot)<rtol) combined2_converging= _TRUE_;
    //printf("I_comb2 = %e, rerr = %e\n",I_comb2,fabs(err/I));
  }


  /* Search for the minimal Laguerre quadrature rule: */
  NLag_max = MIN(N_max,80);
  for (NLag=NL; NLag<=NLag_max; NLag = MIN(NLag_max,NLag+10)){
    /* Evaluate integral: */
    compute_Laguerre(x,w,NLag,0.0,b,c,_TRUE_);
    ILag = 0.0;
    for (i=0; i<NLag; i++){
      (*test)(params_for_function,x[i],&y);
      (*function)(params_for_function,x[i],&y2);
      w[i] *= y2;
      ILag += y*w[i];
    }
    err = I-ILag;
    //fprintf(stderr,"\n Computing Laguerre, N=%d, I=%g and err=%g.\n",NLag,ILag,err);
    if (fabs(err/I)<rtol){
      Laguerre_converging = _TRUE_;
      break;
    }
    Nold = NLag;
    if (Nold == NLag_max) break;
  }

  if (Laguerre_converging == _TRUE_){
    /* We must refine NLag: */
    NL = Nold;
    NR = NLag;
    while ((NR-NL)>1) {
      NLag = (NL+NR)/2;
      /* Evaluate integral: */
      compute_Laguerre(x,w,NLag,0.0,b,c,_TRUE_);
      ILag = 0.0;
      for (i=0; i<NLag; i++){
	(*test)(params_for_function,x[i],&y);
	(*function)(params_for_function,x[i],&y2);
	w[i] *= y2;
	ILag += y*w[i];
      }
      err = I-ILag;
      //fprintf(stderr,"\n NLag=%d, rerr=%g.\n",NLag,fabs(err/I));
      if (fabs(err/Itot)<rtol){
	NR = NLag;
      }
      else{
	NL = NLag;
      }
    }
  }

  /* Choose best method if both works: */
  *N = N_max;
  //Laguerre_converging = _FALSE_;
  if (adapt_converging == _TRUE_) {
    *N = Nadapt;
  }
  if (combined_converging == _TRUE_){
    if(N_comb <= *N){
      *N = N_comb;
      adapt_converging = _FALSE_;
    }
    else{
      combined_converging = _FALSE_;
    }
  }
  if (Laguerre_converging == _TRUE_){
    if (NLag <= *N){
      *N = NLag;
      combined_converging = _FALSE_;
      adapt_converging = _FALSE_;
    }
    else{
      Laguerre_converging = _FALSE_;
    }
  }
  //printf("N_adapt=%d, N_combined=%d at level=%d, Nlag=%d\n",Nadapt,N_comb,level,NLag);
  if (adapt_converging==_TRUE_){
    class_sprintf(method_chosen,"Adaptive Gauss-Kronrod Quadrature");
    /* Gather weights and xvalues from tree: */
    i = Nadapt-1;
    get_leaf_x_and_w(root,&i,x,w,_TRUE_);
  }
  else if (Laguerre_converging==_TRUE_){
    class_sprintf(method_chosen,"Gauss-Laguerre Quadrature");
    /* x and w is already populated in this case. */
  }
  else if (combined_converging == _TRUE_){
    class_sprintf(method_chosen,"Combined Quadrature");
    for(i=0; i<N_comb_leg; i++){
      x[i] = q_leg[i];
      w[i] = w_leg[i];
    }
    i = N_comb_leg;
    get_leaf_x_and_w(root_comb,&i,x,w,_FALSE_);
    //printf("from %d to %d\n",N_comb_leg,i);

    for(i=0; i<N_comb_lag; i++){
      x[N_comb-N_comb_lag+i] = q_lag[i];
      w[N_comb-N_comb_lag+i] = w_lag[i];
    }
  }
  else{
    /* Failed to converge! */
    class_stop(errmsg,
		"get_qsampling fails to obtain a relative tolerance of %g as required using atmost %d points. If the PSD is interpolated from a file, try increasing the resolution and the q-interval (qmin;qmax) if possible, or decrease tol_ncdm/tol_ncdm_bg. As a last resort, increase _QUADRATURE_MAX_/_QUADRATURE_MAX_BG_.",rtol,N_max);
  }
  /* Trim weights to avoid zero weights: */
  for(i=0,zeroskip=0; i<*N; i++){
    for( ;(i<*N)&&(w[i+zeroskip]==0.0); zeroskip++,(*N)--);
    x[i] = x[i+zeroskip];
    w[i] = w[i+zeroskip];
  }


  //printf("Chosen sampling: %s, with %d points.\n",method_chosen,*N);
  //for(i=0; i<*N; i++) printf("(q,w) = (%g,%g)\n",x[i],w[i]);
  /* Deallocate tree: */
  burn_tree(root);
  if(qvec!=NULL){
    burn_tree(root_comb);
    free(wcomb2);
  }

  free(b);
  free(c);

  return _SUCCESS_;
}

int sort_x_and_w(double *x, double *w, double *workx, double *workw, int startidx, int endidx){
  int i,top=endidx,bot=startidx;
  double pivot;
  /* End recursion if only one element left in array: */
  if ((endidx-startidx)<1){
    return _SUCCESS_;
  }
  else{
    /*Copy x and w to workarray: */
    for (i=startidx; i<=endidx; i++){
      workx[i] = x[i];
      workw[i] = w[i];
    }
    pivot = x[endidx];
    //printf("pivot chosen: x[%d] = %g\n",endidx,pivot);
    for (i=startidx; i<endidx; i++){
      if (workx[i]<=pivot){
	//printf("<--%g  ",workx[i]);
	x[bot] = workx[i];
	w[bot++] = workw[i];
      }
      else{
	//printf("  %g-->",workx[i]);
	x[top] = workx[i];
	w[top--] = workw[i];
      }
    }
    //printf("\n top=%d, bot=%d, left=%d, right=%d\n",top,bot,startidx,endidx);
    x[top] = pivot;
    w[top] = workw[endidx];
    /* Recursive call: */
    sort_x_and_w(x,w,workx,workw,startidx,bot-1);
    sort_x_and_w(x,w,workx,workw,top+1,endidx);
    return _SUCCESS_;
  }

}

int get_leaf_x_and_w(qss_node *node, int *ind, double *x, double *w,int isindefinite){
  /* x and w should be exactly 15*root_node->leafchilds, and a leaf count should have
     been performed. Or perhaps I just use the fact that a leaf won't have children.
     Nah, let me use the leaf-count then. */
  int k;
  if (node->leaf_childs==1){
    for(k=0;k<15;k++){
      x[*ind] = node->x[k];
      w[*ind] = node->w[k];
      if (isindefinite == _TRUE_){
	(*ind)--;
      }
      else{
	(*ind)++;
      }
    }
  }
  else{
    /* Do recursive call: */
    get_leaf_x_and_w(node->left,ind,x,w,isindefinite);
    get_leaf_x_and_w(node->right,ind,x,w,isindefinite);
  }
  return _SUCCESS_;
}

int reduce_tree(qss_node *node, int level){
  /* Reduce the tree to a given level. Make all nodes with
     node->leaf_childs==level into leafs.
     If we call reduce_tree(root,1), nothing happens.*/
  if(node->leaf_childs==level){
    burn_tree(node->left);
    burn_tree(node->right);
    node->left = NULL;
    node->right = NULL;
  }
  else if(node->leaf_childs>level){
    /* else try to see if children nodes can be simplified: */
    reduce_tree(node->left,level);
    reduce_tree(node->right,level);
  }
  /* If called on a node which has leaf_childs<level, it does nothing. */
  return _SUCCESS_;
}


int burn_tree(qss_node *node){
  /* Burn node and all subnodes. */
  /* Call burn_branch recursively on children nodes: */
  /* This node and all its subnodes */
  if (node!=NULL){
    if (node->left!=NULL) burn_tree(node->left);
    if (node->right!=NULL) burn_tree(node->right);

    if (node->x!=NULL) free(node->x);
    if (node->w!=NULL) free(node->w);
    free(node);
  }
  return _SUCCESS_;
}

int leaf_count(qss_node *node){
  /* Count the amount of leafs under a given node and write the number in the node. */
  /* We call recursively, until a node is a leaf - then we add the numbers on our
     way back:*/
  if (node->left!=NULL){
    /* This is not a leaf, do recursive call: */
    leaf_count(node->left);
    leaf_count(node->right);
    node->leaf_childs = node->left->leaf_childs + node->right->leaf_childs;
    return _SUCCESS_;
  }
  else{
    /* This is a leaf, by definition leaf_childs = 1: */
    node->leaf_childs = 1;
    return _SUCCESS_;
  }
}

double get_integral(qss_node *node, int level){
  /* Traverse the tree and return the estimate of the integral at a given level.
     level 1 is the best estimate. */
  double IL,IR;
  /* An updated leaf_count is assumed. */
  if (node->leaf_childs<=level){
    return node->I;
  }
  else{
    IL = get_integral(node->left, level);
    IR = get_integral(node->right, level);
    /* Combine the integrals: */
    return (IL+IR);
  }
}



int gk_adapt(
	     qss_node** node,
	     int (*test)(void * params_for_function, double q, double *psi),
	     int (*function)(void * params_for_function, double q, double *f0),
	     void * params_for_function,
	     double tol,
	     int treemode,
	     double a,
	     double b,
	     int isindefinite,
	     ErrorMsg errmsg){
  /* Do adaptive Gauss-Kronrod quadrature, while building the
     recurrence tree. If treemode!=0, store x-values and weights aswell.
     At first call, a and b should be 0 and 1 if isdefinite==_TRUE_. */
  double mid;
  /* Allocate current node: */
  class_alloc(*node,sizeof(qss_node),errmsg);
  if (treemode==0){
    (*node)->x = NULL;
    (*node)->w = NULL;
  }
  else{
    class_alloc((*node)->x,15*sizeof(double),errmsg);
    class_alloc((*node)->w,15*sizeof(double),errmsg);
  }
  (*node)->left = NULL; (*node)->right = NULL;

  gk_quad((*test), (*function), params_for_function, *node, a, b, isindefinite);
  if ((fabs((*node)->err/(*node)->I) < tol)||(tol>=1.0)){
    /* Stop recursion and return. tol>=1.0 in case of I=0 infinite recursion */
    return _SUCCESS_;
  }
  else{
    /* Call gk_adapt recursively on children:*/
    mid = 0.5*(a+b);
    //printf("<-%g,%g,%g,%g",mid,tol,(*node)->err,(*node)->I);
    gk_adapt(&((*node)->left),(*test),(*function), params_for_function, 1.5*tol,
	     treemode, a, mid, isindefinite, errmsg);
    //printf("%g->",mid);
    gk_adapt(&((*node)->right),(*test),(*function), params_for_function, 1.5*tol,
	     treemode, mid, b, isindefinite, errmsg);
    /* Update integral and error in this node and return: */
    /* Actually, it is more convenient just to keep the nodes own estimate of the
       integral for our purposes.
       (*node)->I = (*node)->left->I + (*node)->right->I;
       (*node)->err = sqrt(pow(node->left->err,2)+pow(node->right->err,2));
    */
    return _SUCCESS_;
  }
}

int compute_Hermite(double *x, double *w, int N, int alpha, double *b, double *c){
  int NLag,i;
  double alpha_Lag;

  NLag = N/2;
  /* In case N is uneven, zero the N'th weight:*/
  w[N-1] = 0.0;
  alpha_Lag = (alpha-1.0)/2.0;

  /* Compute the positive roots and weights (up to some simple manipulation): */
  compute_Laguerre(x+NLag,w+NLag,NLag,alpha_Lag,b,c,_FALSE_);

  /* Do manipulations:*/
  for(i=NLag; i<2*NLag; i++){
    x[i] = sqrt(x[i]);
    w[i] *=0.5;
  }

  /* Set the negative roots and weights:*/
  for(i=0; i<NLag; i++){
    x[i] = -x[2*NLag-i-1];
    w[i] = w[2*NLag-i-1];
    if (alpha%2!=0){
      w[i] = -w[i];
    }
  }
  return _SUCCESS_;
}


int compute_Laguerre(double *x, double *w, int N, double alpha, double *b, double *c,int totalweight){
  int i,j,iter,maxiter=10;
  double x0=0.,r1,r2,ratio,d,logprod,logcc;
  double p0,p1,p2,dp0,dp1,dp2;
  double eps=1e-14;
  /* Initialise recursion coefficients: */
  for(i=0; i<N; i++){
    b[i] = alpha + 2.0*i +1.0;
    c[i] = i*(alpha+i);
  }
  logprod = 0.0;
  for(i=1; i<N; i++) logprod +=log(c[i]);
  logcc = lgamma(alpha+1)+logprod;

  /* Loop over roots: */
  for (i=0; i<N; i++){
    /* Estimate root: */
    if (i==0) {
      x0 =(1.0+alpha)*(3.0+0.92*alpha)/( 1.0+2.4*N+1.8*alpha);
    }
    else if (i==1){
      x0 += (15.0+6.25*alpha)/( 1.0+0.9*alpha+2.5*N);
    }
    else{
      r1 = (1.0+2.55*(i-1))/( 1.9*(i-1));
      r2 = 1.26*(i-1)*alpha/(1.0+3.5*(i-1));
      ratio = (r1+r2)/(1.0+0.3*alpha);
      x0 += ratio*(x0-x[i-2]);
    }
    /* Refine root using Newtons method: */
    for(iter=1; iter<=maxiter; iter++){
      /* We need to find p2=L_N(x0), dp2=L'_N(x0) and
	 p1 = L_(N-1)(x0): */
      p1 = 1.0;
      dp1 = 0.0;
      p2 = x0 - alpha - 1.0;
      dp2 = 1.0;
      for (j=1; j<N; j++ ){
	p0 = p1;
	dp0 = dp1;
	p1 = p2;
	dp1 = dp2;
	p2  = (x0-b[j])*p1 - c[j]*p0;
	dp2 = (x0-b[j])*dp1 + p1 - c[j]*dp0;
      }
      /* New guess at root: */
      d = p2/dp2;
      x0 -= d;
      if (fabs(d)<=eps*(fabs(x0)+1.0)) break;
    }
    /* Okay, write root and weight: */
    x[i] = x0;

    if (totalweight == _TRUE_)
      w[i] = exp(x0+logcc-log(dp2*p1));
    else
       w[i] = exp(logcc-log(dp2*p1));
  }

  return _SUCCESS_;

}




int gk_quad(int (*test)(void * params_for_function, double q, double *psi),
	    int (*function)(void * params_for_function, double q, double *f0),
	    void * params_for_function,
	    qss_node* node,
	    double a,
	    double b,
	    int isindefinite){
  const double z_k[15]={-0.991455371120813,
			-0.949107912342759,
			-0.864864423359769,
			-0.741531185599394,
			-0.586087235467691,
			-0.405845151377397,
			-0.207784955007898,
			0.0,
			0.207784955007898,
			0.405845151377397,
			0.586087235467691,
			0.741531185599394,
			0.864864423359769,
			0.949107912342759,
			0.991455371120813};
  const double w_k[15]={0.022935322010529,
			0.063092092629979,
			0.104790010322250,
			0.140653259715525,
			0.169004726639267,
			0.190350578064785,
			0.204432940075298,
			0.209482141084728,
			0.204432940075298,
			0.190350578064785,
			0.169004726639267,
			0.140653259715525,
			0.104790010322250,
			0.063092092629979,
			0.022935322010529};
  const double w_g[7]={0.129484966168870,
		       0.279705391489277,
		       0.381830050505119,
		       0.417959183673469,
		       0.381830050505119,
		       0.279705391489277,
		       0.129484966168870};
  int i,j;
  double x,wg,wk,t,Ik,Ig,y,y2;

  /* 	Loop through abscissas, transform the interval and form the Kronrod
     15 point estimate of the integral.
     Every second time we update the Gauss 7 point quadrature estimate. */

  Ik=0.0;
  Ig=0.0;
  for (i=0;i<15;i++){
    /* Transform z into t in interval between a and b: */
    t = 0.5*(a*(1-z_k[i])+b*(1+z_k[i]));
    /* Modify weight such that it reflects the linear transformation above: */
    wk = 0.5*(b-a)*w_k[i];
    if (isindefinite==_TRUE_){
      /* Transform t into x in interval between 0 and inf: */
      x = 1.0/t-1.0;
      /* Modify weight accordingly: */
      wk = wk/(t*t);
    }
    else{
      x = t;
    }
    (*test)(params_for_function,x,&y);
    (*function)(params_for_function,x,&y2);
    wk *= y2;
    /* Update Kronrod integral: */
    Ik +=wk*y;
    /* If node->x and node->w is allocated, store values: */
    if (node->x!=NULL) node->x[i] = x;
    if (node->w!=NULL) node->w[i] = wk;
    /* If i is uneven, update Gauss integral: */
    if ((i%2)==1){
      j = (i-1)/2;
      /* Transform weight according to linear transformation: */
      wg = 0.5*(b-a)*w_g[j];
      if (isindefinite == _TRUE_){
        /* Transform weight according to non-linear transformation x = 1/t -1: */
        wg = wg/(t*t);
      }
      /* Update integral: */
      Ig +=wg*y*y2;
    }
  }
  node->err = pow(200*fabs(Ik-Ig),1.5);
  node->I = Ik;
  return _SUCCESS_;
}

/**
 * This routine computes the weights and abscissas of a Gauss-Legendre quadrature between -1 and 1
 *
 * @param mu     Input/output: Vector of cos(beta) values
 * @param w8     Input/output: Vector of quadrature weights
 * @param n      Input       : Number of quadrature points
 * @param tol    Input       : tolerance on each mu
 *
 * From Numerical recipes
 **/

int quadrature_gauss_legendre(
			      double *mu,
			      double *w8,
			      int n,
			      double tol,
			      ErrorMsg error_message) {

  int m,j,i,counter;
  double z1,z,pp,p3,p2,p1;

  m=(n+1)/2;
  for (i=1;i<=m;i++) {
    z=cos(_PI_*((double)i-0.25)/((double)n+0.5));
    counter=0;
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
        p3=p2;
        p2=p1;
        p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
      counter++;
      class_test(counter == _MAX_IT_,
		 error_message,
		 "maximum number of iteration reached: increase either _MAX_IT_ or tol\n");
    } while (fabs(z-z1) > tol);
    mu[i-1]=-z;
    mu[n-i]=z;
    w8[i-1]=2.0/((1.0-z*z)*pp*pp);
    w8[n-i]=w8[i-1];

  }
  return _SUCCESS_;
}

int quadrature_in_rectangle(
			    double xl,
			    double xr,
			    double yl,
			    double yr,
			    int *n,
			    double ** x,
			    double ** y,
			    double ** w,
			    ErrorMsg error_message) {


  double xl_tile,xr_tile,yl_tile,yr_tile;
  int N;

  N=24;
  xl_tile = xl;
  xr_tile = xr;
  yl_tile = yl;
  yr_tile = yr;

  *n = N;


  class_alloc(*x,sizeof(double)*N,error_message);
  class_alloc(*y,sizeof(double)*N,error_message);
  class_alloc(*w,sizeof(double)*N,error_message);
  class_call(cubature_order_eleven(xl_tile,
				   xr_tile,
				   yl_tile,
				   yr_tile,
				   *x+0,
				   *y+0,
				   *w+0,
				   error_message),
	     error_message,
	     error_message);


  return _SUCCESS_;
}



int cubature_order_eleven(
			  double xl,
			  double xr,
			  double yl,
			  double yr,
			  double *x,
			  double *y,
			  double *w,
			  ErrorMsg error_message){

  double wi[6]={0.48020763350723814563e-01,
		0.66071329164550595674e-01,
		0.97386777358668164196e-01,
		0.21173634999894860050e+00,
		0.22562606172886338740e+00,
		0.35115871839824543766e+00};
  double xi[6]={0.98263922354085547295e+00,
		0.82577583590296393730e+00,
		0.18858613871864195460e+00,
		0.81252054830481310049e+00,
		0.52532025036454776234e+00,
		0.41658071912022368274e-01};
  double yi[6]={0.69807610454956756478e+00,
		0.93948638281673690721e+00,
		0.95353952820153201585e+00,
		0.31562343291525419599e+00,
		0.71200191307533630655e+00,
		0.42484724884866925062e+00};

  int idx,i;
  double a1,a2,b1,b2;

  a1 = 2./(xr-xl);
  a2 = 2./(yr-yl);
  b1 = 1.-2*xr/(xr-xl);
  b2 = 1.-2*yr/(yr-yl);

  for (i=0,idx=0; i<6; i++,idx++){
    // Upper right corner:
    x[idx] = (xi[i]-b1)/a1;
    y[idx] = (yi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }
  for (i=0,idx=6; i<6; i++,idx++){
    // Upper left corner:
    x[idx] = (-yi[i]-b1)/a1;
    y[idx] = (xi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }
  for (i=0,idx=12; i<6; i++,idx++){
    // Lower left corner:
    x[idx] = (-xi[i]-b1)/a1;
    y[idx] = (-yi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }
  for (i=0,idx=18; i<6; i++,idx++){
    // Lower right corner:
    x[idx] = (yi[i]-b1)/a1;
    y[idx] = (-xi[i]-b2)/a2;
    w[idx] = wi[i]/a1/a2;
  }



  return _SUCCESS_;
}

