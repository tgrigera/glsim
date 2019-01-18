/*
 * timecorr.hh -- Functions for discretized time correlations
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use glsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * glsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#ifndef TIMECORR_HH
#define TIMECORR_HH

#include "fft.hh"

namespace glsim {

/** @{ \ingroup TCorr */

/** 
This routine below computes the time autocorrelation for real
observables. The correlation loop loads [[corr]], the
correlation function is [[c[i]/(N-i)]]. The correlation loop is
$O(N^2)$, and is a direct implementation of the definition. In this
way the elements of [[a]] are read sequentially, so that it turns out
that this loop is slightly more efficient than the alternative using
an outer loop in the time difference (subindex of the [[c]] array)
[Allen].
*/
void correlation_1d_tti_direct(const std::vector<double> &a,
			       std::vector<double> &corr,int nt);

/**
 This function computes the connected correlation. It is computed in
the same way in [[correlation_1d_tti_direct]], except that the average
is substracted to each data point. The average is precomputed in a
separate loop. Although it is possible to compute the average
simultaneously with the nonconnected correlation and then compute the
connected one substracting the square of the average, the present
method introduces less numerical error. Additionally we allow the
caller to request normalizing the correlation to $C(0)=1$ by passing
the boolean [[normalize]].
*/
void correlation_connected_1d_tti_direct(const std::vector<double> &a,
					 std::vector<double> &corr,int nt,
					 bool normalize=false);


/** Compute the correlation at fixed tw by the direct method.
 */
void correlation_1d_tw_direct(const std::vector<double> &a,
			      std::vector<double> &corr,int tw);
/** Compute the connectod correlation at fixed ts by the direct
(O(N^2)) method.
*/
void correlation_connected_1d_tw_direct(const std::vector<double> &a,
					std::vector<double> &corr,int tw,
					bool normalize);

/** Computes self-correlation of a complex series by the direct method */
void correlation_1d_tti_direct(const vcomplex &a,vcomplex &corr,int nt);

/** Computes the connected self-correlation of a complex series by the direct method */
void correlation_connected_1d_tti_direct(const vcomplex &a,
					 vcomplex &corr,int nt,
					 bool normalize=false);

/** Computes self-correlation of scalar real data using FFT (\f$O(N\log N)\f$).

  \param[in] a   Data (time-series)
  \param[out] corr  Correlation function (in the tdata of the FFT object)
  \param[in] nt   Number of correlation points to return (maximum a.size(), suggested 
                  a.size()/2 since points near the end have larger errors)

 */
  void correlation_1d_tti_fft(const std::vector<double> &a,glsim::RealFFT &corr,
			    int nt);

/** Computes connectd self-correlation of scalar real data using FFT (\f$O(N\log N)\f$).

  \param[in] a   Data (time-series)
  \param[out] corr  Correlation function (in the tdata of the FFT object)
  \param[in] nt   Number of correlation points to return (maximum a.size(), suggested 
                  a.size()/2 since points near the end have larger errors)
  \param[in] normalize  Whether to return normalized correlation (i.e. \f$C(0)=1\f$)

 */
void correlation_connected_1d_tti_fft(const std::vector<double> &a,
				      glsim::RealFFT &corr,int nt,
				      bool normalize);


/** Computes self-correlation of scalar complex data using FFT (\f$O(N\log N)\f$).

  \param[in] a   Data (time-series)
  \param[out] corr  Correlation function (in the tdata of the FFT object)
  \param[in] nt   Number of correlation points to return (maximum a.size(), suggested 
                  a.size()/2 since points near the end have larger errors)

 */
void correlation_1d_tti_fft(const vcomplex &a,glsim::ComplexFFT &corr,int nt);

/** Computes connected self-correlation of scalar complex data using FFT (\f$O(N\log N)\f$).

  \param[in] a   Data (time-series)
  \param[out] corr  Correlation function (in the tdata of the FFT object)
  \param[in] nt   Number of correlation points to return (maximum a.size(), suggested 
                  a.size()/2 since points near the end have larger errors)
  \param[in] normalize  Whether to return normalized correlation (i.e. \f$C(0)=1\f$)

 */
void correlation_connected_1d_tti_fft(vcomplex &a,glsim::ComplexFFT &corr,int nt,
				      bool normalize=false);

} /* namespace */

/** @} */

#endif /* TIMECORR_HH */
