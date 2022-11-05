/* specfunc/gsl_sf_hyperg.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_HYPERG22_H__
#define __GSL_SF_HYPERG22_H__

#include <gsl/gsl_sf_result.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

int
hyperg_2F2_series(const double a, const double b, const double c, const double d,
                  const double x, 
                  gsl_sf_result * result
                  );
int
gsl_sf_hyperg_2F2_e(double a, double b, const double c, const double d,
                       const double x,
                       gsl_sf_result * result);

double gsl_sf_hyperg_2F2(double a, double b, double c, double d, double x);

__END_DECLS

#endif /* __GSL_SF_HYPERG_H__ */
