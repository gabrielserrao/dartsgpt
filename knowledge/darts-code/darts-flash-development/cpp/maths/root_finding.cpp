#include <cmath>
#include <vector>
#include <complex>
#include <numeric>
#include <functional>

#include "dartsflash/global/global.hpp"
#include "dartsflash/maths/maths.hpp"
#include "dartsflash/maths/root_finding.hpp"
#include <Eigen/Dense>

RootFinding::RootFinding() {}

int RootFinding::bisection(std::function<double(double)> obj_fun, double& x_, double& a, double& b, double tol_f, double tol_t)
{
   int iter = 1;
   double iter_type = 0;
   (void) iter_type;

   double fa = -1., fb = 1.;
   x = (!std::isnan(x_)) ? x_ : (a+b)*0.5;
   double fx = obj_fun(x);
   if (std::isnan(fx))
   {
      x = NAN;
      return 1;
   }

   if (fx > 0)
   {
      b = x;
      fb = fx;
      fa = (-fb > -1.) ? -1. : -fb;
   }
   else
   {
      a = x;
      fa = fx;
      fb = (-fa < 1.) ? 1. : -fa;
   }

   if (std::fabs(fb) > tol_f)
   {
      // a is the previous value of b and [b, c] always contains the zero.
      double c = a;
      double fc = fa;
      double d = b - a; // the current step
      double e = b - a; // the step at the last step
      (void) e;

      while (true)
      {
         if (fb * fc > 0)  // Check if signs of f(b) and f(c) are opposite
         {
            c = a;
            fc = fa;
            d = b - a;
            e = b - a;
         }
         
         double tol_m = 2 * std::numeric_limits<double>::epsilon() * std::pow(b, 2) + 1e-14;

         // Midpoint.
         double m = (c - b) * 0.5;

         if (std::fabs(fb) < tol_f ||   // objective function has converged - return 0
            (std::abs(m) < tol_m || std::abs(2*m) < tol_t))  // root has converged - discontinuity in objective function, return -1
         {
            x = b;
            return (std::fabs(fb) < tol_f) ? 0 : -1;
         }

         // if (std::abs(e) < tol_m || std::abs(fa) <= std::abs(fb))
         {  // bisection
            d = m;
            e = m;
            iter_type = 1;
         }

         a = b;
         fa = fb;

         if (std::abs(d) > tol_m)
         {
            b = b + d;
         }
         else if (m > 0)
         {
            b = b + tol_m;
         }
         else
         {
            b = b - tol_m;
         }

         fb = obj_fun(b);

         iter++;
         if (std::isnan(fb) || iter > 100)
         {
            x = NAN;
            return 1;
         }
      }
   }

   x = b;
   return 0;
}

int RootFinding::bisection_newton(std::function<double(double)> obj_fun, std::function<double(double)> gradient,
                                  double& x_, double& a, double& b, double tol_f, double tol_t)
{
   (void) gradient;
   int iter = 1;
   double iter_type = 0;
   (void) iter_type;

   double fa = -1., fb = 1.;
   x = (!std::isnan(x_)) ? x_ : (a+b)*0.5;
   double fx = obj_fun(x);
   if (std::isnan(fx))
   {
      x = NAN;
      return 1;
   }

   if (fx > 0)
   {
      b = x;
      fb = fx;
      fa = (-fb > -1.) ? -1. : -fb;
   }
   else
   {
      a = x;
      fa = fx;
      fb = (-fa < 1.) ? 1. : -fa;
   }

   if (std::fabs(fb) > tol_f)
   {
      // a is the previous value of b and [b, c] always contains the zero.
      double c = a;
      double fc = fa;
      double d = b - a; // the current step
      // double e = b - a; // the step at the last step
      // (void) e;

      while (true)
      {
         if (fb * fc > 0)  // Check if signs of f(b) and f(c) are opposite
         {
            c = a;
            fc = fa;
            d = b - a;
            // e = b - a;
         }
         
         double tol_m = 2 * std::numeric_limits<double>::epsilon() * std::pow(b, 2) + 1e-14;

         // Midpoint.
         double m = (c - b) * 0.5;
         double db = -fb/gradient(b);

         if (std::fabs(fb) < tol_f ||   // objective function has converged - return 0
            (std::abs(m) < tol_m || std::abs(2*m) < tol_t))  // root has converged - discontinuity in objective function, return -1
         {
            x = b;
            return (std::fabs(fb) < tol_f) ? 0 : -1;
         }

         // if (std::abs(e) < tol_m || std::abs(fa) <= std::abs(fb))
         std::vector<double> bounds = (b > c) ? std::vector<double>{c, b} : std::vector<double>{b, c};
         if ((db < 0. && b + db > bounds[0]) || (db > 0. && b + db < bounds[1]))
         {
            // Newton step
            d = db;
         }
         else
         {  // bisection step
            d = m;
            // e = m;
            iter_type = 1;
         }
         
         a = b;
         fa = fb;

         if (std::abs(d) > tol_m)
         {
            b = b + d;
         }
         else if (m > 0)
         {
            b = b + tol_m;
         }
         else
         {
            b = b - tol_m;
         }

         fb = obj_fun(b);

         iter++;
         if (std::isnan(fb) || iter > 100)
         {
            x = NAN;
            return 1;
         }
      }
   }

   x = b;
   return 0;
}

int RootFinding::brent(std::function<double(double)> obj_fun, double& x_, double& a, double& b, double tol_f, double tol_t)
{
   // ****************************************************************************
   //
   //  Purpose:
   //
   //    The brent function uses the Brent method to find a root of a given
   //    function f within a given interval [a, b] using an initial estimate t.
   //
   //  Discussion:
   //
   //    The Brent method is an iterative root-finding algorithm that combines the
   //    bisection method,  secant method, and inverse quadratic interpolation to f
   //    ind a root of a given function. It is generally an robust algorithm without
   //    evaluating the derivatives.
   //
   //    The current algorithm is slightly modified based on Moncorge (2022) so that
   //    only one function evaluation is needed to start the algorithm. After the
   //    modification, interval [a, b] should be well-defined to guarantee the root
   //    is bracketed.
   //
   //    Please refer to the original paper of Brent (1971) for other details:
   //    An algorithm with guaranteed convergence for finding a zero of a function.
   //
   //  Parameters:
   //
   //    Input, double a, b, two doubles representing the interval in which to
   //    search for the root.
   //
   //    Input, double t, doubles representing the initial estimate of the
   //    root in the interval [a, b].
   //
   //    Input, double tol_f, tol_t, two doubles representing the tolerances for
   //    the function and the root.
   //
   //    Output, double b, one double representing an approximation of the root of
   //    f within the interval [a, b].
   //
   // ****************************************************************************80

   int iter = 1;
   double iter_type = 0;
   (void) iter_type;

   double fa = -1., fb = 1.;
   x = (!std::isnan(x_)) ? x_ : (a+b)*0.5;
   double fx = obj_fun(x);
   if (std::isnan(fx))
   {
      x = NAN;
      return 1;
   }

   if (fx > 0)
   {
      // Check if the zero is inside [a, b]
      fb = fx;
      fa = (-fb > -1.) ? -1. : -fb;
      b = x;
   }
   else
   {
      // Check if the zero is inside [a, b]
      fa = fx;
      fb = (-fa < 1.) ? 1. : -fa;
      a = x;
   }
   
   if (std::fabs(fb) > tol_f)
   {
      // a is the previous value of b and [b, c] always contains the zero.
      double c = a;
      double fc = fa;
      double d = b - a; // the current step
      double e = b - a; // the step at the last step

      while (true)
      {
         if (fb * fc > 0)  // Check if signs of f(b) and f(c) are opposite
         {
            c = a;
            fc = fa;
            d = b - a;
            e = b - a;
         }

         // Swap to insure f(b) is the smallest value so far.
         if (std::abs(fc) < std::abs(fb))
         {
            a = b;
            fa = fb;
            b = c;
            fb = fc;
            c = a;
            fc = fa;
         }

         double tol_m = 2 * std::numeric_limits<double>::epsilon() * std::pow(b, 2) + 1e-14;

         // Midpoint.
         double m = (c - b) * 0.5;

         if (std::fabs(fb) < tol_f ||   // objective function has converged - return 0
            (std::abs(m) < tol_m || std::abs(2*m) < tol_t))  // root has converged - discontinuity in objective function, return -1
         {
            x = b;
            return (std::fabs(fb) < tol_f) ? 0 : -1;
         }

         if (std::abs(e) < tol_m || std::abs(fa) <= std::abs(fb))
         {  // bisection
            d = m;
            e = m;
            iter_type = 1;
         }
         else
         {
            double pp, qq, rr;   // intermediate variables
            double s = fb / fa;
            if (a == c)
            {  // secant (interpolation)
               pp = 2 * m * s;
               qq = 1 - s;
               iter_type = 2;
            }
            else
            {  // extrapolation
               qq = fa / fc;
               rr = fb / fc;
               pp = s * (2 * m * qq * (qq - rr) - (b - a) * (rr - 1));
               qq = (qq - 1) * (rr - 1) * (s - 1);
               iter_type = 3;
            }

            if (pp > 0)
            {
               qq = -qq;
            }
            else
            {
               pp = -pp;
            }

            s = e;
            e = d;

            if (2 * pp < 3 * m * qq - std::abs(tol_m * qq) && pp < std::abs(0.5 * s * qq))
            {
               // good short step
               d = pp / qq;
            }
            else
            { // bisection
               d = m;
               e = m;
               iter_type = 1;
            }
         }

         a = b;
         fa = fb;

         if (std::abs(d) > tol_m)
         {
            b = b + d;
         }
         else if (m > 0)
         {
            b = b + tol_m;
         }
         else
         {
            b = b - tol_m;
         }

         fb = obj_fun(b);

         iter++;
         if (std::isnan(fb) || iter > 100)
         {
            x = NAN;
            return 1;
         }
      }
   }

   x = b;
   return 0;
}

int RootFinding::brent_newton(std::function<double(double)> obj_fun, std::function<double(double)> gradient, 
                              double& x_, double& a, double& b, double tol_f, double tol_t)
{
   // ****************************************************************************
   //
   //  Purpose:
   //
   //    The brent function uses the Brent method to find a root of a given
   //    function f within a given interval [a, b] using an initial estimate t.
   //
   //  Discussion:
   //
   //    The Brent method is an iterative root-finding algorithm that combines the
   //    bisection method,  secant method, and inverse quadratic interpolation to f
   //    ind a root of a given function. It is generally an robust algorithm without
   //    evaluating the derivatives.
   //
   //    The current algorithm is slightly modified based on Moncorge (2022) so that
   //    only one function evaluation is needed to start the algorithm. After the
   //    modification, interval [a, b] should be well-defined to guarantee the root
   //    is bracketed.
   //
   //    Please refer to the original paper of Brent (1971) for other details:
   //    An algorithm with guaranteed convergence for finding a zero of a function.
   //
   //  Parameters:
   //
   //    Input, double a, b, two doubles representing the interval in which to
   //    search for the root.
   //
   //    Input, double t, doubles representing the initial estimate of the
   //    root in the interval [a, b].
   //
   //    Input, double tol_f, tol_t, two doubles representing the tolerances for
   //    the function and the root.
   //
   //    Output, double b, one double representing an approximation of the root of
   //    f within the interval [a, b].
   //
   // ****************************************************************************
   int iter = 1;
   double iter_type = 0;
   (void) iter_type;

   double fa = -1., fb = 1.;
   double dfa{ NAN }, dfb{ NAN }, d2fb{ NAN };
   x = (!std::isnan(x_)) ? x_ : (a+b)*0.5;
   double fx = obj_fun(x);
   double dfx = gradient(x);
   if (std::isnan(fx))
   {
      x = NAN;
      return 1;
   }

   if (fx > 0)
   {
      // Check if the zero is inside [a, b]
      fb = fx;
      dfb = dfx;
      fa = (-fb > -1.) ? -1. : -fb;
      b = x;
   }
   else
   {
      // Check if the zero is inside [a, b]
      fa = fx;
      dfa = dfx;
      fb = (-fa < 1.) ? 1. : -fa;
      a = x;
   }
   
   if (std::fabs(fb) > tol_f)
   {
      // a is the previous value of b and [b, c] always contains the zero.
      double c = a;
      double fc = fa;
      double dfc = dfa;
      double d = b - a; // the current step
      double e = b - a; // the step at the last step

      while (true)
      {
         if (fb * fc > 0)  // Check if signs of f(b) and f(c) are opposite
         {
            c = a;
            fc = fa;
            dfc = dfa;
            d = b - a;
            e = b - a;
         }

         // Swap to insure f(b) is the smallest value so far.
         if (std::abs(fc) < std::abs(fb))
         {
            a = b;
            fa = fb;
            dfa = dfb;
            b = c;
            fb = fc;
            dfb = dfc;
            c = a;
            fc = fa;
            dfc = dfa;
         }

         double tol_m = 2 * std::numeric_limits<double>::epsilon() * std::pow(b, 2) + 1e-14;

         // Midpoint.
         double m = (c - b) * 0.5;
         double db = -fb / dfb;

         if (std::fabs(fb) < tol_f ||   // objective function has converged - return 0
            (std::abs(m) < tol_m || std::abs(2*m) < tol_t))  // root has converged - discontinuity in objective function, return -1
         {
            x = b;
            return (std::fabs(fb) < tol_f) ? 0 : -1;
         }

         if (iter > 0)
         {
            d2fb = 0.5 * ((dfa-dfb)/(a-b) + (dfc-dfb)/(c-b));
         }

         std::vector<double> bounds = (b > c) ? std::vector<double>{c, b} : std::vector<double>{b, c};
         if (fb * d2fb > 0 && ((db < 0. && b + db > bounds[0]) || (db > 0. && b + db < bounds[1])))
         {  // Newton step
            d = db;
            e = db;
         }
         else if (std::abs(e) < tol_m || std::abs(fa) <= std::abs(fb))
         {  // bisection
            d = m;
            e = m;
            iter_type = 1;
         }
         else
         {
            double pp, qq, rr;   // intermediate variables
            double s = fb / fa;
            if (a == c)
            {  // secant (interpolation)
               pp = 2 * m * s;
               qq = 1 - s;
               iter_type = 2;
            }
            else
            {  // extrapolation
               qq = fa / fc;
               rr = fb / fc;
               pp = s * (2 * m * qq * (qq - rr) - (b - a) * (rr - 1));
               qq = (qq - 1) * (rr - 1) * (s - 1);
               iter_type = 3;
            }

            if (pp > 0)
            {
               qq = -qq;
            }
            else
            {
               pp = -pp;
            }

            s = e;
            e = d;

            if (2 * pp < 3 * m * qq - std::abs(tol_m * qq) && pp < std::abs(0.5 * s * qq))
            {
               // good short step
               d = pp / qq;
            }
            else
            { // bisection
               d = m;
               e = m;
               iter_type = 1;
            }
         }

         a = b;
         fa = fb;
         dfa = dfb;

         if (std::abs(d) > tol_m)
         {
            b = b + d;
         }
         else if (m > 0)
         {
            b = b + tol_m;
         }
         else
         {
            b = b - tol_m;
         }

         fb = obj_fun(b);
         dfb = gradient(b);

         iter++;
         if (std::isnan(fb) || iter > 500)
         {
            x = NAN;
            return 1;
         }
      }
   }

   x = b;
   return 0;
}
