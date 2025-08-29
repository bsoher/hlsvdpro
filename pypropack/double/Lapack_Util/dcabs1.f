*> \brief \b DCABS1_H
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DCABS1_H(Z)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 Z
*       ..
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DCABS1_H computes |Re(.)| + |Im(.)| of a double complex number
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] Z
*> \verbatim
*>          Z is COMPLEX*16
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup double_blas_level1
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION DCABS1_H(Z)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      COMPLEX*16 Z
*     ..
*     ..
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DIMAG
*
      DCABS1_H = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      RETURN
      END