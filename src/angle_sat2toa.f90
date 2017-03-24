SUBROUTINE angle_sat2toa ( earth_curv, ers2_alt, atm_hght, n_arr_pix, sza, vza, relazm )

  USE OMSAO_precision_module,  ONLY: i4, r4, r8
  USE OMSAO_parameters_module, ONLY: &
       r4_missval, deg2rad, rad2deg, pi, min_zenith, max_zenith
  IMPLICIT NONE

  ! Input variables
  ! ===============
  INTEGER (KIND=i4), INTENT (IN) :: n_arr_pix
  REAL    (KIND=r4), INTENT (IN) :: earth_curv, ers2_alt, atm_hght

  ! ==================
  ! Modified variables
  ! ==================
  REAL (KIND=r4), DIMENSION (n_arr_pix), INTENT (INOUT) :: sza, vza, relazm

  ! ===============
  ! Local variables
  ! ===============
  REAL (KIND=r4)                        :: pihalf
  REAL (KIND=r4), DIMENSION (n_arr_pix) :: &
       earth_zen, locazm, szarad, vzarad, razmrad, tmp_sza, tmp_vza


  ! ==================================
  ! Compute GOME viewing angles at TOA
  ! ==================================
  ! (Adopted in part from R.Spurr's geometry module)

  pihalf = REAL(pi/2.0_r8, KIND=r4)


  ! -------------------------------------------------------
  ! Modify a local copy of SZA and VZA so that we don't go
  ! out-of-bounds for the computation of TOA angles.
  ! -------------------------------------------------------
  tmp_sza = sza  ;  tmp_vza = vza
  WHERE ( tmp_sza < min_zenith .OR. tmp_sza > max_zenith )
     tmp_sza = r4_missval
  ENDWHERE
  WHERE ( tmp_vza < min_zenith .OR. tmp_vza > max_zenith )
     tmp_vza = r4_missval
  ENDWHERE

  ! ----------------------------------
  ! First convert everything to Radian
  ! ----------------------------------
  szarad  = sza    * REAL ( deg2rad, KIND=r4 )
  vzarad  = vza    * REAL ( deg2rad, KIND=r4 )
  razmrad = relazm * REAL ( deg2rad, KIND=r4 )
 
  ! ----------------------------------
  ! Step 1: Compute Earth center angle
  ! ----------------------------------    
  earth_zen = ASIN( SIN (vzarad) * (earth_curv+ers2_alt) / (earth_curv+atm_hght) )
  WHERE ( earth_zen <= pihalf )
     earth_zen = earth_zen - vzarad
  ELSEWHERE
     earth_zen = r4_missval
  END WHERE

  ! ------------------------------------------------------------
  ! Step 2: Correction of relative azimuth angle, save as LOCAZM
  ! ------------------------------------------------------------
  WHERE ( earth_zen /= r4_missval )
     locazm = COS(earth_zen) * SIN(szarad) * COS(razmrad) - SIN(earth_zen) * COS(szarad)
  ELSEWHERE
     locazm = r4_missval
  END WHERE

  ! -----------------------------------
  ! Step 3: Correction of zenith angles
  ! -----------------------------------
  WHERE ( earth_zen /= r4_missval )
     szarad = ACOS( COS(szarad) * COS(earth_zen) + SIN(earth_zen) * SIN(szarad) * COS(razmrad) )
     vzarad  = vzarad + earth_zen
  ELSEWHERE
     szarad = r4_missval
     vzarad = r4_missval
  END WHERE
  ! ---------------------------------------------------
  ! Step 4: Final computation of relative azimuth angle
  ! ---------------------------------------------------
  WHERE ( (szarad /= r4_missval) .AND. (SIN(szarad) /= 0.0_r4 ) )
     razmrad = locazm / SIN(szarad)
  ELSEWHERE
     razmrad = r4_missval
  END WHERE

  ! -----------------------------------------------------------------
  ! Make sure the COS of the relative azimuth angle is within [-1,+1]
  ! (i.e., beware of rounding errors), and then take ACOS.
  ! -----------------------------------------------------------------
  WHERE ( (razmrad > 1.0_r4) .AND. (razmrad /= r4_missval) ) 
     razmrad = 1.0_r4
  ENDWHERE
  WHERE ( (razmrad < -1.0_r4) .AND. (razmrad /= r4_missval) ) 
     razmrad = -1.0_r4
  ENDWHERE
  razmrad = ACOS ( razmrad )

  ! ----------------------------------------------
  ! Final step: Convert everything back to Degrees
  ! ----------------------------------------------
  WHERE ( ABS(szarad) <= pihalf )
     sza   = szarad * REAL ( rad2deg, KIND=r4 )
  ELSEWHERE
     sza = r4_missval
  END WHERE
  WHERE ( ABS(vzarad) <= pihalf )
     vza    = vzarad  * REAL ( rad2deg, KIND=r4 )
  ELSEWHERE
     vza  = r4_missval
  ENDWHERE
  WHERE ( ABS(razmrad) <= REAL(pi, KIND=r4) )
     relazm = razmrad * REAL ( rad2deg, KIND=r4 )
  ELSEWHERE
     relazm = r4_missval
  END WHERE

  RETURN
END SUBROUTINE angle_sat2toa
