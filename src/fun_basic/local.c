#include <stdio.h>                          /** for warning in local_divide **/
#include <math.h>
#include <stdlib.h>
/*#include <ieeefp.h> *for finite, on some computers*/

#include "../kernel/coco.h"  /**for DOUBLE**/


/******************************* local_sum ***********************************/

DOUBLE local_sum (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 + val2);
}

/******************************* local_sum_const *****************************/

DOUBLE local_sum_const (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 + par[0]);
}

/******************************* local_sub ***********************************/
/* Attention: the order of the arguments has changed since local_diff!!      */

DOUBLE local_sub (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 - val2);
}

/******************************* local_mult **********************************/

DOUBLE local_mult (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 * val2);
}

/******************************* local_mult_const ****************************/

DOUBLE local_mult_const (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 * par[0]);
}

/******************************* local_div ***********************************/
/* formerly local_divide                                                     */

DOUBLE local_div (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (val2 != 0.0)
        return (val1 / val2);
    else {
        static int firsttime = 1;

        if  (val1 != 0.0)
            if  (firsttime) {
                fprintf (stderr, "\n\nlocal_div: dividing nonzero number by zero -- not warning again!\n\n");
                firsttime = 0;
            }

        return (0.0);
    }
}

/******************************* local_div_const *****************************/

DOUBLE local_div_const (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (par[0] != 0.0)
        return (val1 / par[0]);
    else {
        if  (val1 != 0.0)
            fprintf (stderr, "\nlocal_div_const divide by zero!  ");
        return (0.0);
    }
}

/******************************* local_modulo_const **************************/

DOUBLE local_modulo_const (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    DOUBLE retval = val1;

    while (retval >= par[0])
        retval -= par[0];

    while (retval < 0.0)                 /**!!new!!**/
        retval += par[0];

    return retval;
}

/******************************* local_cyclic ********************************/
/* needed after taking difference between two cyclic variables ...           */

DOUBLE local_cyclic (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    DOUBLE retval = val1;

    while (retval >= 0.5 * par[0])
        retval -= par[0];

    while (retval < -0.5 * par[0])
        retval += par[0];

    return retval;
}

/******************************* local_average *******************************/
/* Parameter 0 <= par[0](aux1) <= 1. !!                                      */

DOUBLE local_average (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 * par[0] + val2 * (1.0 - par[0]));
}

/******************************* local_copy **********************************/

DOUBLE local_copy (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1);
}

/******************************* local_persist *******************************/

DOUBLE local_persist (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  ((val1 == 0) && (val2 > 0.0))
        return (0.0);
    else
        return (1.0);
}

/******************************* local_threshold *****************************/

DOUBLE local_threshold (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (fabs(val1) < par[0])
        return (0.0);
    else
        return (val1);
}

/******************************* local_isbigger ******************************/

DOUBLE local_isbigger (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (val1 > val2)
        return (1.0);
    else
        return (0.0);
}

/******************************* local_biggerof ******************************/

DOUBLE local_biggerof (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 > val2 ? val1 : val2);
}

/******************************* local_biggerfabsof **************************/

DOUBLE local_biggerfabsof (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (fabs(val1) > fabs(val2) ? val1 : val2);
}




/******************************* local_const *********************************/
/* Replaces local_zero now with par[0] = 0.                                  */

DOUBLE local_const (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (par[0]);
}

/******************************* local_rectify *******************************/

DOUBLE local_rectify (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (par[1] == 1.0)
        return (val1 < par[0] ? par[0] : val1);
    else /**-1.0**/
        return (val1 > par[0] ? par[0] : val1);
}

/******************************* local_invert ********************************/

DOUBLE local_invert (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (-val1);
}

/******************************* local_abs ***********************************/

DOUBLE local_abs (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (fabs (val1));
}

/******************************* local_sign **********************************/

DOUBLE local_sign (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 >= 0.0 ? 1.0 : -1.0);
}

/******************************* local_round *********************************/

DOUBLE local_round (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (par[0] == 0.0)
        fprintf (stderr, "\nlocal_round wants a non-zero parameter, e.g. 0.5");

    return (val1 > par[0] ? 1.0 : 0.0);
}

/******************************* local_sqrt **********************************/

DOUBLE local_sqrt (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (sqrt (val1));
}

/******************************* local_power *********************************/

DOUBLE local_power (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (pow (fabs (val1), par[0]));
}

/******************************* local_exp ***********************************/

DOUBLE local_exp (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (par[0] == -1)
        return (exp (-val1));
    else
        return (exp (val1));
}

/******************************* local_log_pos *******************************/

DOUBLE local_log_pos (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

  if  (val1 >= 0)
    return (log (val1 + par[0]));
  else
    return (- log (-val1 + par[0]));
}

/******************************* local_cos ***********************************/

DOUBLE local_cos (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (cos (val1));
}

/******************************* local_sin ***********************************/

DOUBLE local_sin (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (sin (val1));
}

/******************************* local_tan ***********************************/

DOUBLE local_tan (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (tan (val1));
}

/******************************* local_atan **********************************/

DOUBLE local_atan (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (atan (val1));
}

/******************************* local_atan_to_2pi ***************************/
/* val1=im (y-axis), val2=re (x-axis)                                        */

DOUBLE local_atan_to_2pi (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    /**y-axis**/
    if  (val2 == 0.0) {
        if  (val1 >= 0.0)
            return 0.5 * M_PI;
        else
            return 1.5 * M_PI;
    }

    /**1st quadrant and pos x-axis**/
    if  ((val1 >= 0.0) && (val2 > 0.0))
        return atan (val1 / val2);

    /**3rd quadrant and neg x-axis**/
    if  ((val1 <= 0.0) && (val2 < 0.0))
        return M_PI + atan (val1 / val2);

    /**2nd quadrant**/
    if  ((val1 > 0.0) && (val2 < 0.0))
        return M_PI + atan (val1 / val2);

    /**4th quadrant**/
    if  ((val1 < 0.0) && (val2 > 0.0))
        return 2.0 * M_PI + atan (val1 / val2);

    fprintf (stderr, "\nlocal_atan_to_2pi should not come to this point!\n");

    return 0;
}

/******************************* local_lin_01 ********************************/

DOUBLE local_lin_01 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 < 0.0 ? 0.0 : (val1 > 1.0 ? 1.0 : val1));
}




/******************************* local_rand **********************************/
/*  formerly: -1 + 2 * drand48()                                             */

DOUBLE local_rand (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (par[0] + (par[1] - par[0]) * drand48());
}

/******************************* local_rand_pos_neg **************************/

DOUBLE local_rand_pos_neg (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (drand48() < 0.5 ? -1.0 : 1.0);
}

/******************************* local_as_rand *******************************/
/* par[0] scales input                                                       */
/* par[1] scales output (like sat)                                           */

DOUBLE local_as_rand (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (drand48() < (par[0] * val1) ? par[1] : 0.0);
}

/******************************* local_circ_gauss ****************************/
/* par[0]=mean, par[1]=sigma. Gaussian with x-axis squashed on circle 0..2PI.*/
/* f(x) = x / (PI - abs (x))                                                 */
/* g(x) = exp (-0.5 * x*x / (sigma*sigma))                                   */
/* plot g(f(x))                                                              */
/* Used by single_circ_gauss.                                                */

DOUBLE local_circ_gauss (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    DOUBLE diff = fabs (val1 - par[0]);

    if  (diff >= M_PI)
        diff = 2.0 * M_PI - diff;

    if  ((diff < 0.0) || (diff > M_PI))
        fprintf (stderr, "\nlocal_circ_gauss: diff out of boundary!\n");

    DOUBLE f = M_PI - diff;

    if  (f == 0.0)
        return (DOUBLE)0;

    f = diff / f;

    return (exp (-0.5 * f * f / (par[1] * par[1])));
}

/******************************* local_rand_gibbs ****************************/
/* Parameter par[0](kurt) is sparseness.                                     */
/* Parameter par[1](sat)  is return value.                                   */

DOUBLE local_rand_gibbs (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (drand48() < (1.0 / (1.0 + par[0]))
            ? (drand48() < 0.5 ? par[1] : -par[1])
            : 0.0);
}

/******************************* local_rand_gibbs_01 *************************/
/* Parameter par[0](kurt) is sparseness.                                     */
/* Parameter par[1](sat)  is return value.                                   */

DOUBLE local_rand_gibbs_01 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (drand48() < (1.0 / (1.0 + par[0])) ? par[1] : 0.0);
}

/******************************* local_rand_exp ******************************/
/* Produces exponentially distributed values >= 0. Parameter par[0] = mean.  */

DOUBLE local_rand_exp (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    /**Exponentialfunktion 1/b * exp(-x/b) von [0..X]**/
    /**Gilt: 1/b int_0^X exp(-x/b) dx = -exp(-x/b) |_0^X = 1 - exp(-X/b) =y **/
    /**Umkehrfunktion davon: X = -b * log(1-y)**/

    return -par[0] * log (1.0 - drand48 ());
}

/******************************* local_tanh **********************************/

DOUBLE local_tanh (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (tanh (val1));
}

/******************************* local_bp_tanh *******************************/
/* Derivative of tanh(R) = 1-tanh(R)^2. To be applied to output vals of tanh.*/

DOUBLE local_bp_tanh (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return ((1.0 - val1 * val1) * val2);
}




/******************************* local_sparse ********************************/
/* returns 1 - d/du S(u) where p(u) = exp (-S(u)); see Diss. P.55            */
/* Parameter par[0](kurt) for sparseness.                                    */

DOUBLE local_sparse (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 - par[0] * 2.0 * val1 / (1.0 + val1 * val1));
}

/******************************* local_sparse_diff ***************************/
/* Derivative of local_sparse (with respect to original inputs).             */
/* par[0](kurt) should be < 0.5, else derivative at 0 <= 0!                  */

DOUBLE local_sparse_diff (DOUBLE *par, DOUBLE val1, DOUBLE val2) {
    DOUBLE term = 1.0 / (1.0 + val1 * val1);

    return (1.0 - 2.0 * par[0] * term
                + 4.0 * par[0] * val1 * val1 * term * term);
}

/******************************* local_sparse_01 *****************************/
/* Parameter par[0](kurt) for sparseness.                                    */

DOUBLE local_sparse_01 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return ((val1 > 0.0) ? val1 - par[0] * 2.0 * val1 / (1.0 + val1 * val1)
                         : 0.0);
}


/******************************* local_zhang *********************************/
/* K.Zhang, J.Neurosci. (1996) Repr. ... head direction cell ensemble.       */
/* Accordingly: bias=0.5  b=10  beta=0.8  a=6.34.  I suggest: bias=0  a=0.2  */

DOUBLE local_zhang (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double bias = par[0];   /**shift ~ threshold**/
    double b    = par[1];   /**scales the input**/
    double beta = par[2];   /**smaller bends log-like -- larger makes linear**/
    double a    = par[3];   /**scales the output**/
    double exponent, inner, out;
    int linearizing = 0;

    exponent = b * (val1 + bias);

    if  (exponent > 707.0)
        linearizing = 1;

    if  (linearizing)
        exponent = 707.0;

    inner = log (1.0 + exp (exponent));

    out = a * pow (inner, beta);

    if  (linearizing)
        fprintf (stderr, "\nlocal_zhang: too large input=%f -- cut output=%f  ", val1, out);

    return (out);
}
/* gnuplot
 bias = 0.0 ; b = 10.0 ; beta = 0.8 ; a = 0.2
 inner (x) = log (1.0 + exp (b * (x + bias)))
 f(x) = a * exp (log (inner(x)) * beta)
*/

/******************************* local_zhang_scale ***************************/
/* Like local_zhang, but par[3] (output scale) is replaced by 2nd arg, val2. */
/* HACKED: NOT SCALE BUT BIAS (THRESHOLD) IS REPLACED!                       */

DOUBLE local_zhang_scale (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double bias = val2;     /**shift ~ threshold**/
    double b    = par[1];   /**scales the input**/
    double beta = par[2];   /**smaller bends log-like -- larger makes linear**/
    double a    = par[3];   /**scales the output**/
    double exponent, inner, out;
    int linearizing = 0;

    exponent = b * (val1 + bias);

    if  (exponent > 707.0)
        linearizing = 1;

    if  (linearizing)
        exponent = 707.0;

    inner = log (1.0 + exp (exponent));

    out = a * pow (inner, beta);

    if  (linearizing)
        fprintf (stderr, "\nlocal_zhang: too large input=%f -- cut output=%f  ", val1, out);

    return (out);
}


/******************************* local_zhang2 ********************************/
/* K.Zhang, J.Neurosci. (1996) Repr. ... head direction cell ensemble.       */
/* No parameter beta! ~like other function if beta=0.35.                     */

DOUBLE local_zhang2 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {
                           /**if param is increased ...                     **/
    double bias = par[0];  /**shifts the curve to the left                  **/  /**shifts all weights down**/
    double b    = par[1];  /**bends, like increasing temp in fermi-function **/  /**scales weight differences down (doesn't change mean which is det by bias**/
    double a    = par[2];  /**scales the curve up.                          **/  /**          "                                 "                           **/

    double inner = log (1.0 + exp (b * (val1 + bias)));

    return (a * log (1.0 + inner));
}


/******************************* local_zhang_inv *****************************/

DOUBLE local_zhang_inv (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double bias = par[0];
    double b    = par[1];
    double beta = par[2];
    double a    = par[3];

    double inner = pow ((val1 / a), (1.0 / beta));

    return ((log (exp (inner) - 1.0) - b * bias) / b);
}
/*maple> solve (a * (ln (1 + exp (b * (x + c))))^beta = sigma, x);
                                      /__1_\
c=bias;                        /sigma\\beta/
sigma is                ln(exp(|-----|      ) - 1) - b c
previous                _______\__a__/__________________
output                                  b
*/

/******************************* local_zhang2_inv ****************************/

DOUBLE local_zhang2_inv (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double bias = par[0];
    double b    = par[1];
    double a    = par[2];

    return ((log (exp (exp (val1 / a) - 1.0) - 1.0) - b * bias) / b);
}
/*maple> solve (a * ln (1 + ln(1 + exp(b * (x + bias)))) = s, x);    
s is                   ln(exp(exp(s/a) - 1) - 1) - b bias
previous               ----------------------------------
output                                  b
*/

/******************************* local_zhang_diff ****************************/

DOUBLE local_zhang_diff (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double bias = par[0];
    double b    = par[1];
    double beta = par[2];
    double a    = par[3];

    double inner = 1.0 + exp (b * (val1 + bias));

    return ( (a * pow (log (inner), beta) * beta * b * exp (b * (val1 + bias)))
           / (inner * log (inner)));
}
/*maple> diff (a * (ln (1 + exp (b * (x + c))))^beta, x);
                                       beta
               a ln(1 + exp(b (x + c)))     beta b exp(b (x + c))
               --------------------------------------------------
                   (1 + exp(b (x + c))) ln(1 + exp(b (x + c)))
*/

/******************************* local_zhang2_diff ***************************/

DOUBLE local_zhang2_diff (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double bias = par[0];
    double b    = par[1];
    double a    = par[2];

    double inner = exp (b * (val1 + bias));

    return ( (a * b * inner)
           / ((1.0 + inner) * (1.0 + log (1.0 + inner))));
}
/*maple> diff (a * ln (1 + ln(1 + exp(b * (x + bias)))), x);
                             a b exp(b (x + bias))
            -------------------------------------------------------
            (1 + exp(b (x + bias))) (1 + ln(1 + exp(b (x + bias))))
*/


/******************************* local_mean **********************************/
/* Mean-field version of ternary (-1, 0, 1) Gibbs.                           */
/* par[0](kurt) for sparseness.                                              */
/* par[1](beta) for inverse temperature.                                     */
/* par[2](sat)  for max amplitude.                                           */
/* f(x) = p2 * (exp (p1*x) - exp (-p1*x)) / (exp (p1*x) + p0 + exp (-p1*x))  */

DOUBLE local_mean (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h     = exp ( par[1] * val1);
    double exp_min_h = exp (-par[1] * val1);

    return (par[2] * (exp_h - exp_min_h) / (exp_h + par[0] + exp_min_h));
}


/******************************* local_bp_mean *******************************/
/* Does not exist because derivative cannot be expressed by the function.    */

/******************************* local_mean_diff *****************************/
/* Derivative of local_mean (with respect to original inputs).               */
/* par[0](kurt) for sparseness.                                              */
/* par[1](beta) must be 1.0!                                                 */

DOUBLE local_mean_diff (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h     = exp ( par[1] * val1);
    double exp_min_h = exp (-par[1] * val1);
    double nenner    = exp_h + par[0] + exp_min_h;

    return (  (exp_h + exp_min_h) / nenner
            - (exp_h - exp_min_h) * (exp_h - exp_min_h) / (nenner * nenner));
    /**
    gnuplot> n = 10.0; f(x) = (exp(x) - exp(-x)) / (exp(x) + n + exp(-x))
    gnuplot>           df(x) = (exp(x) + exp(- x)) / (exp(x) + n + exp(- x))
                       - ((exp(x) - exp(- x))**2) / (exp(x) + n + exp(- x))**2
    **/
}


/******************************* local_mean_01 *******************************/
/* mean-field version of binary (0, 1) Gibbs.                                */
/* par[0](kurt) sparseness.                                                  */
/* par[1](beta) inverse Temperature.                                         */
/* par[2](aux2) resting act (degeneracy of spontaneous ON-state)             */
/* par[3](bias) negative threshold, comparable to local_zhang (new!).        */
/* par[4](sat)  max amplitude.                                               */

DOUBLE local_mean_01 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h = exp (par[1] * (val1 + par[3]));

    if  (finite (exp_h))
        return (par[4] * (exp_h + par[2]) / (exp_h + par[2] + par[0]));
    else
        return (par[4]);
}
/*
gnuplot
p0=1;p1=1;p2=0;p3=0;p4=1
e(x) = exp(p1 * (x + p3))
f(x) = p4 * (e(x) + p2) / (e(x) + p2 + p0)
set xrange [-5:5]
plot f(x) notitle
*/


/******************************* local_mean_scale ****************************/
/* like local_mean_01 but beta (par[1]) is replaced by 2nd arg, val2.        */

DOUBLE local_mean_scale (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h = exp (val2 * (val1 + par[3]));

    if  (finite (exp_h))
        return (par[4] * (exp_h + par[2]) / (exp_h + par[2] + par[0]));
    else
        return (1.0);
}

/******************************* local_mean_scale_diff **************************/
/* like local_mean_01_diff but beta (par[1]) is replaced by 2nd arg, val2.      */

DOUBLE local_mean_scale_diff (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h = exp (val2 * (val1 + par[3]));
    double nenner = exp_h + par[2] + par[0];

    return (par[4] * val2 * exp_h / nenner
          - par[4] * (exp_h + par[2]) * val2 * exp_h / (nenner * nenner));
}


/******************************* local_bp_mean_01 ****************************/
/* Derivative of local_mean_01, to be applied to output values of mean_01 !! */

DOUBLE local_bp_mean_01 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return (val1 - val1 * val1);
}

/******************************* local_mean_01_inv ***************************/
/* Inverse.                                                                  */

DOUBLE local_mean_01_inv (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    return ((log (-(par[4] * par[2] - val1 * par[2] - val1 * par[0]) / (par[4] - val1)) - par[1] * par[3]) / par[1]);
}

/******************************* local_mean_01_diff **************************/
/* Derivative of local_mean_01 (with respect to original inputs).            */

DOUBLE local_mean_01_diff (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h = exp (par[1] * (val1 + par[3]));
    double nenner = exp_h + par[2] + par[0];

    return (par[4] * par[1] * exp_h / nenner
          - par[4] * (exp_h + par[2]) * par[1] * exp_h / (nenner * nenner));
}

/******************************* local_mean_01_d0 ****************************/
/* Derivative of local_mean_01 w.r.t. parameter kurt.                        */

DOUBLE local_mean_01_d0 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h  = exp (par[1] * (val1 + par[3]));
    double nenner = exp_h + par[2] + par[0];

    return (- par[4] * (exp_h + par[2]) / (nenner * nenner));
}

/******************************* local_mean_01_d1 ****************************/
/* Derivative of local_mean_01 w.r.t. parameter beta.                        */

DOUBLE local_mean_01_d1 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h  = exp (par[1] * (val1 + par[3]));
    double nenner = exp_h + par[2] + par[0];

    return ( par[4] * (val1 + par[3]) * exp_h / nenner
           - (par[4] * (exp_h + par[2]) * (val1 + par[3]) * exp_h)
           / (nenner * nenner));
}

/******************************* local_mean_01_d2 ****************************/
/* Derivative of local_mean_01 w.r.t. parameter rest_act.                    */

DOUBLE local_mean_01_d2 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h  = exp (par[1] * (val1 + par[3]));
    double nenner = exp_h + par[2] + par[0];

    return (par[4] / nenner - par[4] * (exp_h + par[2]) / (nenner * nenner));
}

/******************************* local_mean_01_d3 ****************************/
/* Derivative of local_mean_01 w.r.t. parameter bias.                        */

DOUBLE local_mean_01_d3 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h  = exp (par[1] * (val1 + par[3]));
    double nenner = exp_h + par[2] + par[0];

    return ( par[4] * par[1] * exp_h / nenner
           - par[4] * (exp_h + par[2]) * par[1] * exp_h / (nenner * nenner));
}

/******************************* local_mean_01_d4 ****************************/
/* Derivative of local_mean_01 w.r.t. parameter sat.                         */

DOUBLE local_mean_01_d4 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h  = exp (par[1] * (val1 + par[3]));
    double nenner = exp_h + par[2] + par[0];

    return ((exp_h + par[2]) / nenner);
}


/******************************* local_gibbs *********************************/
/* Stochastic ternary (-1, 0, 1) transfer function.                          */
/* par[0](kurt) for sparseness.                                              */
/* par[1](beta) is inverse temperature.                                      */
/* par[2](sat)  for return value.                                            */

DOUBLE local_gibbs (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h     = exp ( par[1] * val1);
    double exp_min_h = exp (-par[1] * val1);
    double part_sum  = exp_h + par[0] + exp_min_h;
    double rd = drand48 ();

    if  (rd < (exp_min_h / part_sum)) {
        return (-par[2]);
    } else {
        if  (rd > (1.0 - exp_h / part_sum))
                                   /**=((exp_min_h+x->kurt_t[ari])/part_sum)**/
            return (par[2]);
        else
            return (0.0);
    }
}

/******************************* local_gibbs_01 ******************************/
/* Stochastic binary (0,1) transfer function.                                */
/* par[0](kurt) for sparseness.                                              */
/* par[1](beta) for inverse temparature.                                     */
/* par[2](aux2) resting activity (degeneracy of spontaneous ON-state!)       */
/* par[3](bias) negative threshold, comparable to local_zhang (new!).        */
/* par[4](sat)  for return value.                                            */

DOUBLE local_gibbs_01 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h = exp (par[1] * (val1 + par[3]));

    if  (drand48 () < (exp_h + par[2]) / (exp_h + par[2] + par[0]))
        return (par[4]);
    else
        return (0.0);
}

/******************************* local_gibbs_01_2 ****************************/
/* Stochastic binary (0,1) transfer function.                                */
/* Sparseness fixed at 1.0 (only one other state).                           */
/* par[0](beta) for inverse temparature.                                     */
/* par[1](sat)  for return value.                                            */

DOUBLE local_gibbs_01_2 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    double exp_h     = exp (par[0] * val1);
    double part_sum  = exp_h + 1.0;
    double rd = drand48 ();

    if  (rd < exp_h / part_sum)           /**(rd > (1.0 - exp_h / part_sum))**/
        return (par[1]);
    else
        return (0.0);
}



/******************************* local_odds_01 *******************************/
/* Stochastic binary (0,1) update function for positive inputs.              */
/* See "The Helmholtz Machine", p.900, for motivation (Dayan, Hinton, Neal). */
/* par[0](sat)  for return value.      No temperature, no kurtosis (yet).    */

DOUBLE local_odds_01 (DOUBLE *par, DOUBLE val1, DOUBLE val2) {

    if  (drand48 () < 1.0 - 1.0 / (1.0 + val1))
        return (par[0]);
    else
        return (0.0);
}
