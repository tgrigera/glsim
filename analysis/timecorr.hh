/*
 * timecorr.hh -- Functions for discretized time correlations
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, with the additional requirements of
 * attribution and nonmisrepresentation. You can use either version 3, or
 * (at your option) any later version.
 * 
 * Additional terms under GNU GPL version 3 section 7:
 * 
 * When you redistribute this software, you are required to preserve its
 * author attributions. If you distribute verbatim copies, you must not
 * alter the AUTHORS file or attributions inserted in the source files,
 * and you must not change the software's name. If you distribute a
 * modified copy, then you must give clear notice that your work is
 * different from but based on glsim. You must distribute it under a
 * different name, but include a prominent notice specifying that "(your
 * package) is based on glsim version x.x", and provide a pointer to the
 * glsim distribution.
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
