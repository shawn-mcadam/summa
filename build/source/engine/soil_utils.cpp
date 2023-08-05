#include <cmath> // std::abs
#include <limits> // std::numeric_limits<double>::epsilon()
#include <LUST/lust.hpp> // lust::pow(double, double)

// constants
static double dx=-1e-12; // finite difference increment

// slightly more convenient to use this than to write std::pow(x,2)
static inline double square(double x){ return x*x; }


extern "C" double hydCond_psi(double psi, double k_sat, double alpha, double n, double m);
double hydCond_psi(double psi, double k_sat, double alpha, double n, double m){
  auto x   = psi*alpha;
  auto xn1 = lust::pow(x, n-1.0);
  if(psi < 0)
    return k_sat*square(1.0 - xn1*lust::pow(1.0 + x*xn1, -m))/lust::pow(1.0 + x*xn1, m/2.0);
  else
    return k_sat;
}


extern "C" double hydCond_liq(double volFracLiq, double k_sat, double theta_res, double theta_sat, double m);
double hydCond_liq(double volFracLiq, double k_sat, double theta_res, double theta_sat, double m){
  if(volFracLiq < theta_sat){
    auto theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res); // effective soil moisture
    return k_sat*lust::pow(theta_e, 0.5) * square(1.0 - lust::pow(1.0 - lust::pow(theta_e, 1.0/m), m));
  } else return k_sat;
}


extern "C" double volFracLiq(double psi, double alpha, double theta_res, double theta_sat, double n, double m);
double volFracLiq(double psi, double alpha, double theta_res, double theta_sat, double n, double m){
  if(psi < 0.0)
    return theta_res + (theta_sat - theta_res)*lust::pow(1.0 + lust::pow(alpha*psi, n), -m);
  else
    return theta_sat;
}


extern "C" double matricHead(double theta, double alpha, double theta_res, double theta_sat, double n, double m);
double matricHead(double theta, double alpha, double theta_res, double theta_sat, double n, double m){
  // we want the effective saturation to be at least epsilon to avoid division by 0
  auto epsilon = std::numeric_limits<double>::epsilon();
  auto effSat  = std::max(epsilon, (theta - theta_res) / (theta_sat - theta_res));

  // effSat is necessarily positive
  return (effSat < 1.0) ? lust::pow(lust::pow(effSat,-1.0/m) - 1.0, 1.0/n) / alpha : 0.0;
}


/* note that psi*alpha > 0 because fortran would return nan otherwise */
extern "C" double dTheta_dPsi(double psi, double alpha, double theta_res, double theta_sat, double n, double m);
double dTheta_dPsi(double psi, double alpha, double theta_res, double theta_sat, double n, double m){
  auto epsilon = std::numeric_limits<double>::epsilon();
  if(psi > 0.0) return epsilon;

  auto x = psi*alpha;
  auto xn1 = lust::pow(x,n-1.0);
  auto dTdP = -(theta_sat-theta_res)*m*n*alpha*xn1/lust::pow(1.0 + x*xn1, m+1.0);
  return (std::abs(dTdP) < epsilon) ? epsilon : dTdP;
}


extern "C" double dPsi_dTheta(double volFracLiq, double alpha, double theta_res, double theta_sat, double n, double m);
double dPsi_dTheta(double volFracLiq, double alpha, double theta_res, double theta_sat, double n, double m){
  // check if less than saturation
  if(volFracLiq >= theta_sat) return 0.0;

  // compute effective water content TODO does max(0.001,...) really make any sense here? Should that be epsilon? or dx?
  auto theta_e = std::max(0.001,(volFracLiq - theta_res) / (theta_sat - theta_res)); // effective soil moisture
  // compute the 1st double and derivative
  auto pow_theta_e1 = lust::pow(theta_e, -1.0/m - 1);
  auto y1 = pow_theta_e1*theta_e - 1;
  auto d1 = -pow_theta_e1 / (m*(theta_sat - theta_res));
  // compute the 2nd double and derivative
  auto powy1 = lust::pow(y1, 1.0/n - 1);
  auto y2 = powy1*y1;
  auto d2 = powy1/n;
  return d1*d2/alpha;
}


extern "C" double dPsi_dTheta2(double volFracLiq, double alpha, double theta_res, double theta_sat, double n, double m, bool lTangent);
double dPsi_dTheta2(double volFracLiq, double alpha, double theta_res, double theta_sat, double n, double m, bool lTangent){
  // return 0 if the volumetric liquid water content exceeds porosity
  if(volFracLiq >= theta_sat) return 0.0;

  if(lTangent){ // ***** compute derivative using the analytical formula
    // theta_e is the effective saturation (effective soil moisture)
    auto theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res);

    // get the first double and derivative
    auto mthroot_theta_e = lust::pow(theta_e,1.0/m);
    auto y1 = -1.0 / (m*(theta_sat - theta_res)*mthroot_theta_e*theta_e);
    auto d1 = (m + 1) / (square(m*(theta_sat-theta_res))*mthroot_theta_e*theta_e*theta_e);
    // get the second double and derivative
    auto xx = 1.0 / mthroot_theta_e - 1; auto nthroot_xx2 = lust::pow(xx, 1/n-2);
    auto y2 = nthroot_xx2*xx/n;
    auto d2 = (n - 1)*nthroot_xx2 / ((theta_sat - theta_res)*m*n*n* mthroot_theta_e*theta_e);
    return (d1*y2 + y1*d2)/alpha;
  }else{ // ***** approximate derivative with forward difference
    auto f0 = dPsi_dTheta(volFracLiq,   alpha,theta_res,theta_sat,n,m);
    auto f1 = dPsi_dTheta(volFracLiq+dx,alpha,theta_res,theta_sat,n,m);
    return (f1 - f0)/dx;
  }
}


extern "C" double dHydCond_dPsi(double psi, double k_sat, double alpha, double n, double m, bool lTangent);
double dHydCond_dPsi(double psi, double k_sat, double alpha, double n, double m, bool lTangent){
  // derivative is zero if saturated
  if(psi >= 0.0) return 0.0;

  // ***** compute analytical derivatives
  if(lTangent){
    // compute the derivative for the numerator
    auto x = psi*alpha;
    auto xn2 = lust::pow(x, n - 2);
    auto xn1 = xn2*x;
    auto xn0 = xn2*x*x;

    auto f_x1 =  xn1;
    auto f_x2 =  lust::pow(1.0 + xn0, -m);
    auto d_x1 =  alpha*(n - 1)*xn2;
    auto d_x2 = -alpha*m*n* xn1 *lust::pow(1 + xn0, -m - 1.0);
    auto f_nm =  square(1.0 - f_x1*f_x2);
    auto d_nm = -2*(d_x1*f_x2 + f_x1*d_x2) * (1 - f_x1*f_x2);
    // compute the derivative for the denominator
    auto f_dm = lust::pow(1.0 + xn0, m/2.0);
    auto d_dm = 0.5*alpha*m*n*xn1*lust::pow(1 + xn0, m/2.0 - 1.0);
    // and combine
    return k_sat*(d_nm*f_dm - d_dm*f_nm) / (f_dm*f_dm);
  }else{
    // ***** approximate derivative with forward difference
    auto f0 = hydCond_psi(psi,   k_sat,alpha,n,m);
    auto f1 = hydCond_psi(psi+dx,k_sat,alpha,n,m);
    return (f1 - f0)/dx;
  }
}


extern "C" double dHydCond_dLiq(double volFracLiq, double k_sat, double theta_res, double theta_sat, double m, bool lTangent);
double dHydCond_dLiq(double volFracLiq, double k_sat, double theta_res, double theta_sat, double m, bool lTangent){
  // derivative is zero if super-saturated
  if(volFracLiq >= theta_sat) return 0.0;

  // ***** compute analytical derivatives
  if(lTangent){
    // compute the effective saturation
    auto theta_e = (volFracLiq - theta_res) / (theta_sat - theta_res);
    // compute the double and derivative of the first fuction
    auto sqrt_theta_e = lust::pow(theta_e, 0.5);
    auto f1 = k_sat*sqrt_theta_e;
    auto d1 = k_sat*0.5/(sqrt_theta_e*(theta_sat - theta_res));
    // compute the double and derivative of the second double
    // (first part)
    auto mthroot_theta_e1 = lust::pow(theta_e, 1/m - 1);
    auto x1 = 1 - mthroot_theta_e1*theta_e;
    auto p1 = -mthroot_theta_e1/(m*(theta_sat - theta_res));  // differentiate (1.d - theta_e**(1.d/m)
    // (second part)
    auto x1m1 = lust::pow(x1, m - 1);
    auto x2 = x1*x1m1;
    auto p2 =  m*x1m1;
    // (final)
    auto f2 = square(1 - x2);
    auto p3 = -2*(1 - x2);
    // (combine)
    auto d2 = p1*p2*p3;
    // pull it all together
    return d1*f2 + d2*f1;
  }else{
    // ***** compute numerical derivatives
    auto f0 = hydCond_liq(volFracLiq,   k_sat,theta_res,theta_sat,m);
    auto f1 = hydCond_liq(volFracLiq+dx,k_sat,theta_res,theta_sat,m);
    return (f1 - f0)/dx;
  }
}


extern "C" double dTheta_dTk(double Tk, double theta_res, double theta_sat, double alpha, double n, double m);
double dTheta_dTk(double Tk, double theta_res, double theta_sat, double alpha, double n, double m){
  /* constants. Ideally, these should be imported from source/dshare/multiconst.f90, but that is a tricky dependency to work out... */
  auto LH_fus  = 333700.0;               // latent heat of fusion   (m2 s-2)
  auto gravity = 9.80616;                // acceleration of gravity (m s-2)
  auto Tfreeze = 273.16;                 // temperature at freezing (K)
  auto kappa = LH_fus/(gravity*Tfreeze); // (m K-1)

  // differentiate the freezing curve w.r.t. temperature -- making use of the chain rule
  auto xtemp = alpha*kappa*(Tk-Tfreeze); // (-)
  auto xtempn1 = lust::pow(xtemp,n-1);
  return -(theta_sat - theta_res)*alpha*kappa*m*n*xtempn1*lust::pow(1.0 + xtemp*xtempn1, -m - 1.0);
}

