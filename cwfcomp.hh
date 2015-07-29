using namespace std;

// Direct integration of the Coulomb equation
// ------------------------------------------
// One uses the Burlisch-Stoer-Henrici method, where one integrates on different meshes
// with the Henrici method, and then uses the Richardson method to get the final result by extrapolation.
// Numerical Recipes, Chap. 16.4 .

class ODE_integration
{
public:
  ODE_integration (const complex<double> &l_1,
		   const complex<double> &two_eta_1)
    : l (l_1), ll_plus_one (l_1*(l_1+1.0)), two_eta (two_eta_1)
    {
      for (int n = 0 ; n < 8 ; n++)
	for (int i = 0 ; i < n ; i++)
	{
	  interpolation_term_tab[n][i] = 1.0;

	  for (int j = 0 ; j < n ; j++)
	    if (i != j)
	      interpolation_term_tab[n][i] *= (i+1.0)/(i-j);
	}

      for (unsigned int k = 0 ; k < 8 ; k++) m_tab[k] = 2*(k+1);
      for (unsigned int k = 0 ; k < 8 ; k++) one_over_m_tab[k] = 1.0/static_cast<double> (m_tab[k]);
    }

  void operator() (const complex<double> &r0,const complex<double> &u0,const complex<double> &du0,
		   const complex<double> &r,complex<double> &u,complex<double> &du) const;

  private:
  complex<double> extrapolation_in_zero (const unsigned int n,const complex<double> *f) const;
  complex<double> F_r_u (const complex<double> &r,const complex<double> &u) const;
  void integration_Henrici (const unsigned int m,const complex<double> &h,
			    const complex<double> &r0,const complex<double> &u0,const complex<double> &du0,
			    const complex<double> &r,complex<double> &u,complex<double> &du) const;
  
  const complex<double> l,ll_plus_one;  // angular momentum,l(l+1).
  const complex<double> two_eta;        // 2.eta, with eta the Sommerfeld parameter.

  unsigned int m_tab[8];                                 // integers used in the extrapolation method.
  double one_over_m_tab[8],interpolation_term_tab[8][8]; // doubles used in the extrapolation method.
};



// Extrapolation in h=0 of a table of function values h close to h=0
// -----------------------------------------------------------------
//
// Variables:
// ----------
// n : number of points of the function f near h=0.
// T : table containing the points f[h(0)]...f[h(n-1)] close to h=0.
// f_in_zero : extrapolated value of the points f[h(0)]...f[h(n-1)] in h=0.

 
complex<double> ODE_integration::extrapolation_in_zero (const unsigned int n,const complex<double> T[]) const
{  
  complex<double> f_in_zero = 0.0;

  for (unsigned int i = 0 ; i < n ; i++)
    f_in_zero += interpolation_term_tab[n][i]*T[i];

  return f_in_zero;
}



// Calculation of F(z,u(z)) in u''(z) = F(z,u(z))  
// ----------------------------------------------
//
// F(z,u(z))=(l(l+1)/(z^2) + 2.eta/z - 1).u(z), 
//
// Variables:
// ----------
// z : parameter of the wave function.
// u : discretized wave function in z.
// one_over_z : 1.0/z

complex<double> ODE_integration::F_r_u (const complex<double> &z,const complex<double> &u) const
{
  if (l == 0) return (two_eta/z - 1.0)*u;

  const complex<double> one_over_z = 1.0/z;

  return ((ll_plus_one*one_over_z + two_eta)*one_over_z - 1.0)*u;
}



// Integration with discretization of u''(r)=F(r,u(r)) with the Henrici method.
// ----------------------------------------------------------------------------
//
// See Numerical Recipes for the method.
//
// Initials conditions : r0,u(r0),du/dr(r0).
// Obtained functions : r,u(r),du/dr(r).
//
// Variables:
// ----------
// m : number of intervals between r0 and r
// h : integration step (r-r0)/m .
// r0,u0,du0 : r0,u(r0),du/dr(r0).
// r,u,du : r,u(r),du/dr(r).
// h_square : h*h
// half_h = 0.5*h
// delta : value used in the Henrici method.

void ODE_integration::integration_Henrici (const unsigned int m,const complex<double> &h,
					   const complex<double> &r0,const complex<double> &u0,const complex<double> &du0,
					   const complex<double> &r,complex<double> &u,complex<double> &du) const
{
  const complex<double> h_square = h*h,half_h = 0.5*h;

  complex<double> delta = h*(du0 + half_h*F_r_u (r0,u0));
  u = u0 + delta;

  for (unsigned int i = 1 ; i < m ; i++)
  {
    delta += h_square*F_r_u (r0 + i*h,u);
    u += delta;
  }
  
  du = delta/h + half_h*F_r_u (r,u);
}










// Integration of u''(r) = F(r,u(r)) with the Bulirsch-Stoer-Henrici method.
// -------------------------------------------------------------------------
//
// Initials conditions : r0,u0=u(r0),du0=du/dr(r0)
// Obtained functions : r,u=u(r),du=du/dr(r)
//
// See Numerical Recipes for the method.
//
// Variables:
// ---------- 
// r0,u0,du0 : r0,u(r0),du/dr(r0).
// r,u,du : r,u(r),du/dr(r).
// H,r_debut,r_end,u_debut,du_debut : length of an integration interval, debut and end of the integration interval, u(r_debut), u'(r_debut).
//                                    H is equal to r-r0 at the beginning and is divided by 2 each time the extrapolation fails with 16 sub-intervals 
//                                    between r_debut and r_end. If H = [r-r0]/N, with N=2,4,8,16,...,
//                                    the integration intervals are [r_debut = r0:r_end = r0+H], ..., [r_debut = r0+(N-1).H,r_end = r].
// u_end,du_end,k : tables of u(r_end) and u'(r_end) values calculated with the Henrici method with 2,4,6,...,2(k+1) sub-intervals between r_debut and r_end,
//                  with 0 <= k <= 7.
// H_over_m_tab : H/m for m=2,4,6,...,16.
// inf_norm_half_H : |H/2|oo. It is used to know if r = r_debut up to numerical accuracy, as one has |r-r_debut|oo <= |H/2|oo for this case only.
// u_extrapolated,u_extrapolated_next : values of u extrapolated from the points of the table u_end with k and k->k+1 points, k >= 2.
// test : test to know if the method worked, i.e., |u_extrapolated/u_extrapolated_next - 1|oo < precision.
// du_extrapolated_next : value of u'(r_end) extrapolated from k points of the table du_end, k >= 3.
// r_debut_plus_H : r_debut+H. r_debut+H at the end of integration is not necessarily r because of numerical cancellations.
//                             In this case, r_end must be put equal to r.

void ODE_integration::operator() (const complex<double> &r0,const complex<double> &u0,const complex<double> &du0,
				  const complex<double> &r,complex<double> &u,complex<double> &du) const
{
  if (r == r0) {u = u0; du = du0; return;}

  complex<double> r_debut = r0,u_debut = u0,du_debut = du0,H = r-r0,u_end[8],du_end[8],u_extrapolated_next,du_extrapolated_next;
  double test = 1.0;

  while (test > precision)
  {
    complex<double> H_over_m_tab[8];
    for (unsigned int k = 0 ; k < 8 ; k++) H_over_m_tab[k] = H*one_over_m_tab[k];
    const double inf_norm_half_H = inf_norm (H_over_m_tab[0]);

    while (inf_norm (r_debut - r) > inf_norm_half_H)
    {
      const complex<double> r_debut_plus_H = r_debut + H, r_end = (inf_norm (r - r_debut_plus_H) > inf_norm_half_H) ? (r_debut_plus_H) : (r);

      integration_Henrici (2,H_over_m_tab[0],r_debut,u_debut,du_debut,r_end,u_end[0],du_end[0]);
      integration_Henrici (4,H_over_m_tab[1],r_debut,u_debut,du_debut,r_end,u_end[1],du_end[1]);
      complex<double> u_extrapolated = extrapolation_in_zero (2,u_end); 

      unsigned int k = 2; 
      do
      {
	integration_Henrici (m_tab[k],H_over_m_tab[k],r_debut,u_debut,du_debut,r_end,u_end[k],du_end[k]);
	u_extrapolated_next = extrapolation_in_zero (++k,u_end);
	test = inf_norm (u_extrapolated/u_extrapolated_next - 1.0);
	u_extrapolated = u_extrapolated_next;
      }
      while ((test > precision) && (k < 7));

      r_debut += H;
      u_debut = u_extrapolated_next;
      du_debut = du_extrapolated_next = extrapolation_in_zero (k,du_end);
    }

    H *= 0.5;
    r_debut = r0;
    u_debut = u0;
    du_debut = du0;
  }

  u = u_extrapolated_next;
  du = du_extrapolated_next;
}











// Class to calculate the Coulomb wave functions
// ---------------------------------------------

class Coulomb_wave_functions
{
public:

  // Constructor.
  // ------------
  // Constants are defined in the constructor, 
  // plus a pointer to class ODE_integration, ODE_ptr, to integrate numerically the regular Coulomb wave function.
  // 
  // Variables:
  // ----------
  // is_it_normalized_c : true if one wants normalized functions, i.e. the standard normalization,
  //                      false if one wants F -> F/C(l,eta) and H+/H-/G -> H+/H-/G.C(l,eta), to avoid overflows for |eta| >> 1 and |z| small.
  // l_c : orbital angular momentum.
  // eta_c : Sommerfeld parameter.

  Coulomb_wave_functions (const bool is_it_normalized_c,const complex<double> &l_c,const complex<double> &eta_c)
    : l (l_c),
      is_it_normalized (is_it_normalized_c),
      eta (eta_c), 
      neg_int_omega_one ((rint (real (l_c + complex<double> (-imag (eta_c),real (eta_c)))) == l_c + complex<double> (-imag (eta_c),real (eta_c))) && 
			 (rint (real (1 + l_c + complex<double> (-imag (eta_c),real (eta_c)))) <= 0.0)),
      neg_int_omega_minus_one ((rint (real (l_c - complex<double> (-imag (eta_c),real (eta_c)))) == l_c - complex<double> (-imag (eta_c),real (eta_c))) && 
			       (rint (real (1 + l_c - complex<double> (-imag (eta_c),real (eta_c)))) <= 0.0)),
      sigma_l (sigma_l_calc (l_c,eta_c)),
      log_Cl_eta (log_Cl_eta_calc (l_c,eta_c)),
      Cl_eta (exp (log_Cl_eta_calc (l_c,eta_c))),
      exp_I_chi (exp_I_omega_chi_calc (1,l_c,eta_c)),
      exp_minus_I_chi (exp_I_omega_chi_calc (-1,l_c,eta_c)),
      one_over_sin_chi (1.0/sin_chi_calc (l_c,eta_c)),
      log_cut_constant_CFa_plus (log_cut_constant_CFa_calc (is_it_normalized_c,1,l_c,eta_c)),
      log_cut_constant_CFa_minus (log_cut_constant_CFa_calc (is_it_normalized_c,-1,l_c,eta_c)),
      cut_constant_CFa_plus (exp (log_cut_constant_CFa_calc (is_it_normalized_c,1,l_c,eta_c))),
      cut_constant_CFa_minus (exp (log_cut_constant_CFa_calc (is_it_normalized_c,-1,l_c,eta_c))),
      log_cut_constant_CFb_plus (log_cut_constant_CFb_calc (is_it_normalized_c,1,l_c,eta_c)),
      log_cut_constant_CFb_minus (log_cut_constant_CFb_calc (is_it_normalized_c,-1,l_c,eta_c)),
      log_cut_constant_AS_plus (log_cut_constant_AS_calc (1,l_c,eta_c)),
      log_cut_constant_AS_minus (log_cut_constant_AS_calc (-1,l_c,eta_c)),
      cut_constant_CFb_plus (exp (log_cut_constant_CFb_calc (is_it_normalized_c,1,l_c,eta_c))),
      cut_constant_CFb_minus (exp (log_cut_constant_CFb_calc (is_it_normalized_c,-1,l_c,eta_c))),
      log_sym_constant_arg_neg ((is_it_normalized_c) ? (-M_PI*(eta_c+(l_c+1)*complex<double> (0.0,1.0))) : (-M_PI*(l_c+1)*complex<double> (0.0,1.0))),
      log_sym_constant_arg_pos ((is_it_normalized_c) ? (-M_PI*(eta_c-(l_c+1)*complex<double> (0.0,1.0))) : (M_PI*(l_c+1)*complex<double> (0.0,1.0))),
      sym_constant_arg_neg ((is_it_normalized_c) ? (exp (-M_PI*(eta_c+(l_c+1)*complex<double> (0.0,1.0)))) : (exp (-M_PI*(l_c+1)*complex<double> (0.0,1.0)))),
      sym_constant_arg_pos ((is_it_normalized_c) ? (exp (-M_PI*(eta_c-(l_c+1)*complex<double> (0.0,1.0)))) : (exp (M_PI*(l_c+1)*complex<double> (0.0,1.0)))), 
      turning_point (max (1.0,abs (eta_c) + sqrt (abs (l_c*(l_c+1.0)) + abs (eta_c*eta_c)))),
      is_H_dir_int_naive (false),cwf_real_ptr (0),cwf_real_eta_plus_ptr (0),cwf_real_eta_minus_ptr (0),cwf_real_l_plus_ptr (0),cwf_real_l_minus_ptr (0),
      cwf_minus_eta_ptr (0),cwf_lp_ptr (0),prec_first_order_expansion (0.1*sqrt_precision)
  {
    ODE_ptr = new class ODE_integration (l,2.0*eta);

    debut = 0.0;

    if (real (l) >= 0.0)
    {
      F_debut = 0.0;
      dF_debut = (l == 0) ? ((is_it_normalized) ? (Cl_eta) : (1.0)) : (0.0);
    }
  }

  ~Coulomb_wave_functions (void)
  {
    delete cwf_real_ptr;

    delete cwf_real_l_plus_ptr;
    delete cwf_real_l_minus_ptr;

    delete cwf_real_eta_plus_ptr;
    delete cwf_real_eta_minus_ptr;

    delete cwf_minus_eta_ptr;
    delete cwf_lp_ptr;

    delete ODE_ptr;
  }

  void F_dF_init (const complex<double> &z,const complex<double> &F,const complex<double> &dF);

  void F_dF (const complex<double> &z,complex<double> &F,complex<double> &dF); 
  void G_dG (const complex<double> &z,complex<double> &G,complex<double> &dG);
  void H_dH (const int omega,const complex<double> &z,complex<double> &H,complex<double> &dH);
  void H_dH_scaled (const int omega,const complex<double> &z,complex<double> &H,complex<double> &dH);

  const complex<double> l,eta; // Angular momentum and Sommerfeld parameter.

  const bool is_it_normalized;
  // true if F(z) ~ C(l,eta).z^{l+1} in 0, false if F(z) ~ z^{l+1} in 0.

private:

  void asymptotic_series (const int omega,const complex<double> &one_over_z,complex<double> sum[],complex<double> dsum[],bool &is_it_successful);

  complex<double> continued_fraction_f (const complex<double> &z,const int omega);
  complex<double> continued_fraction_h (const complex<double> &z,const int omega);

  void F_dF_power_series (const complex<double> &z,complex<double> &F,complex<double> &dF);

  void asymptotic_expansion_F_dF (const complex<double> &z,complex<double> &F,complex<double> &dF,bool &is_it_successful);
  void asymptotic_expansion_H_dH_scaled (const int omega,const complex<double> &one_over_z,
					 complex<double> &H_scaled,complex<double> &dH_scaled,bool &is_it_successful);

  void F_dF_direct_integration (const complex<double> &z,complex<double> &F,complex<double> &dF,bool &is_it_successful);
  void H_dH_direct_integration (const int omega,const complex<double> &z,complex<double> &H,complex<double> &dH);

  void partial_derivatives (const bool is_it_regular,const bool is_it_eta,const double x,double &d_chi_Bx,double &d_chi_dBx);
  void first_order_expansions (const bool is_it_regular,const complex<double> &z,complex<double> &B,complex<double> &dB);
  void H_dH_from_first_order_expansions (const int omega,const complex<double> &z,complex<double> &H,complex<double> &dH);

  void H_dH_with_F_dF_and_CF (const int omega,const complex<double> &z,complex<double> &H,complex<double> &dH);
  void H_dH_with_expansion (const int omega,const complex<double> &z,complex<double> &H,complex<double> &dH,bool &is_it_successful);

  void F_dF_with_symmetry_relations (const complex<double> &z,complex<double> &F,complex<double> &dF);

  const bool neg_int_omega_one,neg_int_omega_minus_one;
  // neg_int_omega_one : true if 1+l+i.eta is negative integer, false if not.
  // neg_int_omega_minus_one : true if 1+l-i.eta is negative integer, false if not.

  const complex<double> log_Cl_eta,Cl_eta,sigma_l; // log[C(l,eta)], C(l,eta), sigma(l,eta)

  const complex<double> cut_constant_CFa_plus,cut_constant_CFa_minus,cut_constant_CFb_plus,cut_constant_CFb_minus;
  const complex<double> log_cut_constant_CFa_plus,log_cut_constant_CFa_minus,log_cut_constant_CFb_plus,log_cut_constant_CFb_minus;
  const complex<double> log_cut_constant_AS_plus,log_cut_constant_AS_minus;
  // cut constants and their logs for continued fractions (CFa and CFb) and asymptotic series (AS).
  // plus,minus is for omega = 1 or -1.
  // See functions log_cut_constant_AS_calc, log_cut_constant_CFa_calc and log_cut_constant_CFb_calc.

  const complex<double> exp_I_chi,exp_minus_I_chi,one_over_sin_chi; 
  // exp (i.chi), exp (-i.chi), 1/sin (chi) with chi = sigma(l,eta) - sigma(-l-1,eta) - (l+1/2).Pi .
  // They are used to calculate H+/H- from F(l,eta,z) and F(-l-1,eta,z) if |Im[l]| >= 1 and |z| <= 1 .

  const complex<double> sym_constant_arg_neg,sym_constant_arg_pos,log_sym_constant_arg_neg,log_sym_constant_arg_pos;
  // Multiplicative constants and their logs used in the following reflection formulas : 
  // F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta-i.l)] if arg (z) > 0 and is_it_normalized is true, so sym_constant_arg_pos = -exp[-Pi.(eta-i.l)],
  // F(l,eta,z) = -F(l,-eta,-z).exp[-Pi.(eta+i.l)] if arg (z) <= 0 and is_it_normalized is true, so sym_constant_arg_neg = -exp[-Pi.(eta-i.l)],
  // F(l,eta,z) = -F(l,-eta,-z).exp[i.Pi.l)] if arg (z) > 0 and is_it_normalized is false, so sym_constant_arg_pos = -exp[i.Pi.l)],
  // F(l,eta,z) = -F(l,-eta,-z).exp[-i.Pi.l)] if arg (z) <= 0 and is_it_normalized is false, so sym_constant_arg_neg = -exp[-i.Pi.l)].

  const double turning_point,prec_first_order_expansion; // turning_point : max (1,||eta| + sqrt[|l(l+1)| + |eta|^2]|).
                                                         // prec_first_order_expansion : 0.1*sqrt_precision. It is the precision used for first_order_expansions.

  bool is_H_dir_int_naive; // true if one integrates H+/H- forward without considering |H+/H-|, false if not. It is false except in continued_fraction_h.

  complex<double> debut,F_debut,dF_debut;
  // Coulomb wave functions and derivative at z = debut.
  // It is used to integrate the Coulomb wave function faster, 
  // as debut is usually close to the argument of the Coulomb wave function so that the integration is quicker and more stable.

  class ODE_integration *ODE_ptr;  // pointer to class ODE_integration to integrate numerically the Coulomb equation.

  class Coulomb_wave_functions *cwf_real_ptr,*cwf_real_l_plus_ptr,*cwf_real_l_minus_ptr,*cwf_real_eta_plus_ptr,*cwf_real_eta_minus_ptr;
  // pointers to classes Coulomb_wave_functions of parameters (l_r,eta_r) (one has eta_r = Re[eta], eta_i = Im[eta], l_r = Re[l] and l_i = Im[l]), 
  // (l_r +/- epsilon[l],eta_r) and (l_r,eta_r +/- epsilon[eta]).
  // They are first put to zero and allocated in the program when they are needed.
  // They are used for the first order expansion method when |l_i| << 1, |eta_i| << 1 and |Im[z]| << Re[z] with Re[z] > 0.

  class Coulomb_wave_functions *cwf_minus_eta_ptr,*cwf_lp_ptr;
  // pointers to classes Coulomb_wave_functions of parameters (l,-eta), (lp = -l-1,eta) and (l_r +/- precision,eta).
  // They are first put to zero and allocated in the program when they are needed.
  // They are used for symmetry relations : F(l,eta,z) \propto F(l,-eta,-z) and h[omega](l,eta,z) = -h[omega](l,-eta,-z)
  // and to calculate H+/H- from F(l,eta,z) and F(lp,eta,z) if |Im[l]| >= 1 and |z| <= 1.
};


