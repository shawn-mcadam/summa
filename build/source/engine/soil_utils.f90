! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module soil_utils_module

! data types
use, intrinsic :: iso_c_binding
USE nrtype

USE multiconst,only: gravity, & ! acceleration of gravity       (m s-2)
                     Tfreeze, & ! temperature at freezing    (K)
                     LH_fus,  & ! latent heat of fusion      (J kg-1, or m2 s-2)
                     R_wv       ! gas constant for water vapor  (J kg-1 K-1; [J = Pa m3])

! privacy
implicit none
private

! routines to make public
public::iceImpede
public::dIceImpede_dTemp
public::hydCond_psi
public::hydCond_liq
public::hydCondMP_liq
public::dHydCond_dPsi
public::dHydCond_dLiq
public::volFracLiq
public::matricHead
public::dTheta_dPsi
public::dPsi_dTheta
public::dPsi_dTheta2
public::RH_soilair
public::dTheta_dTk
public::crit_soilT
public::liquidHead
public::gammp

! constant parameters
real(rkind),parameter     :: valueMissing=-9999._rkind    ! missing value parameter

interface
 ! ******************************************************************************************************************************
 ! public function hydCond_psi: compute the hydraulic conductivity as a function of matric head (m s-1)
 ! ******************************************************************************************************************************
 function hydCond_psi(psi,k_sat,alpha,n,m) bind(C, name="hydCond_psi")
 ! computes hydraulic conductivity given psi and soil hydraulic parameters k_sat, alpha, n, and m
 import c_double
 implicit none
 real(c_double),value :: psi         ! soil water suction (m)
 real(c_double),value :: k_sat       ! saturated hydraulic conductivity (m s-1)
 real(c_double),value :: alpha       ! scaling parameter (m-1)
 real(c_double),value :: n           ! vGn "n" parameter (-)
 real(c_double),value :: m           ! vGn "m" parameter (-)
 real(c_double)       :: hydCond_psi ! hydraulic conductivity (m s-1)
 end function hydCond_psi


 ! ******************************************************************************************************************************
 ! public function hydCond_liq: compute the hydraulic conductivity as a function of volumetric liquid water content (m s-1)
 ! ******************************************************************************************************************************
 function hydCond_liq(volFracLiq,k_sat,theta_res,theta_sat,m) bind(C, name="hydCond_liq")
 ! computes hydraulic conductivity given volFracLiq and soil hydraulic parameters k_sat, theta_sat, theta_res, and m
 import c_double
 implicit none
 real(c_double),value :: volFracLiq  ! volumetric liquid water content (-)
 real(c_double),value :: k_sat       ! saturated hydraulic conductivity (m s-1)
 real(c_double),value :: theta_res   ! residual volumetric liquid water content (-)
 real(c_double),value :: theta_sat   ! soil porosity (-)
 real(c_double),value :: m           ! vGn "m" parameter (-)
 real(c_double)       :: hydCond_liq ! hydraulic conductivity (m s-1)
 end function hydCond_liq


 ! ******************************************************************************************************************************
 ! public function volFracLiq: compute the volumetric liquid water content (-)
 ! ******************************************************************************************************************************
 function volFracLiq(psi,alpha,theta_res,theta_sat,n,m) bind(C, name="volFracLiq")
 ! computes the volumetric liquid water content given psi and soil hydraulic parameters theta_res, theta_sat, alpha, n, and m
 import c_double
 implicit none
 real(c_double),value :: psi         ! soil water suction (m)
 real(c_double),value :: alpha       ! scaling parameter (m-1)
 real(c_double),value :: theta_res   ! residual volumetric water content (-)
 real(c_double),value :: theta_sat   ! porosity (-)
 real(c_double),value :: n           ! vGn "n" parameter (-)
 real(c_double),value :: m           ! vGn "m" parameter (-)
 real(c_double)       :: volFracLiq  ! volumetric liquid water content (-)
 end function volFracLiq


 ! ******************************************************************************************************************************
 ! public function matricHead: compute the matric head (m) based on the volumetric liquid water content
 ! ******************************************************************************************************************************
 function matricHead(theta,alpha,theta_res,theta_sat,n,m) bind(C, name="matricHead")
 ! computes the volumetric liquid water content given psi and soil hydraulic parameters theta_res, theta_sat, alpha, n, and m
 import c_double
 implicit none
 real(c_double),value :: theta       ! volumetric liquid water content (-)
 real(c_double),value :: alpha       ! scaling parameter (m-1)
 real(c_double),value :: theta_res   ! residual volumetric water content (-)
 real(c_double),value :: theta_sat   ! porosity (-)
 real(c_double),value :: n           ! vGn "n" parameter (-)
 real(c_double),value :: m           ! vGn "m" parameter (-)
 real(c_double)       :: matricHead  ! matric head (m)
 end function matricHead


 ! ******************************************************************************************************************************
 ! public function dTheta_dPsi: compute the derivative of the soil water characteristic (m-1)
 ! ******************************************************************************************************************************
 function dTheta_dPsi(psi, alpha, theta_res, theta_sat, n, m) bind(C, name="dTheta_dPsi")
 import c_double
 implicit none
 real(c_double),value :: psi         ! soil water suction (m)
 real(c_double),value :: alpha       ! scaling parameter (m-1)
 real(c_double),value :: theta_res   ! residual volumetric water content (-)
 real(c_double),value :: theta_sat   ! porosity (-)
 real(c_double),value :: n           ! vGn "n" parameter (-)
 real(c_double),value :: m           ! vGn "m" parameter (-)
 real(c_double)       :: dTheta_dPsi ! derivative of the soil water characteristic (m-1)
 end function


 ! ******************************************************************************************************************************
 ! public function dPsi_dTheta: compute the derivative of the soil water characteristic (m-1)
 ! ******************************************************************************************************************************
 function dPsi_dTheta(volFracLiq,alpha,theta_res,theta_sat,n,m) bind(C, name="dPsi_dTheta")
 import c_double
 implicit none
 real(c_double),value :: volFracLiq  ! volumetric liquid water content (-)
 real(c_double),value :: alpha       ! scaling parameter (m-1)
 real(c_double),value :: theta_res   ! residual volumetric water content (-)
 real(c_double),value :: theta_sat   ! porosity (-)
 real(c_double),value :: n           ! vGn "n" parameter (-)
 real(c_double),value :: m           ! vGn "m" parameter (-)
 real(c_double)       :: dPsi_dTheta ! derivative of the soil water characteristic (m)
 end function dPsi_dTheta


 ! ******************************************************************************************************************************
 ! public function dPsi_dTheta2: compute the derivative of dPsi_dTheta (m-1)
 ! ******************************************************************************************************************************
 function dPsi_dTheta2(volFracLiq,alpha,theta_res,theta_sat,n,m,lTangent) bind(C, name="dPsi_dTheta2")
 import c_double, c_bool
 implicit none
 real(c_double),value :: volFracLiq   ! volumetric liquid water content (-)
 real(c_double),value :: alpha        ! scaling parameter (m-1)
 real(c_double),value :: theta_res    ! residual volumetric water content (-)
 real(c_double),value :: theta_sat    ! porosity (-)
 real(c_double),value :: n            ! vGn "n" parameter (-)
 real(c_double),value :: m            ! vGn "m" parameter (-)
 logical(c_bool),value :: lTangent     ! method used to compute derivative (.true. = analytical)
 real(c_double)       :: dPsi_dTheta2 ! derivative of the soil water characteristic (m)
 end function dPsi_dTheta2


 ! ******************************************************************************************************************************
 ! public function dHydCond_dPsi: compute the derivative in hydraulic conductivity w.r.t. matric head (s-1)
 ! ******************************************************************************************************************************
 function dHydCond_dPsi(psi,k_sat,alpha,n,m,lTangent) bind(C, name="dHydCond_dPsi")
 ! computes the derivative in hydraulic conductivity w.r.t matric head,
 !  given psi and soil hydraulic parameters k_sat, alpha, n, and m
 import c_double, c_bool
 implicit none
 real(c_double),value :: psi            ! soil water suction (m)
 real(c_double),value :: k_sat          ! saturated hydraulic conductivity (m s-1)
 real(c_double),value :: alpha          ! scaling parameter (m-1)
 real(c_double),value :: n              ! vGn "n" parameter (-)
 real(c_double),value :: m              ! vGn "m" parameter (-)
 logical(c_bool),value :: lTangent       ! method used to compute derivative (.true. = analytical)
 real(c_double)       :: dHydCond_dPsi  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 end function dHydCond_dPsi


 ! ******************************************************************************************************************************
 ! public function dHydCond_dLiq: compute the derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 ! ******************************************************************************************************************************
 ! computes the derivative in hydraulic conductivity w.r.t the volumetric fraction of liquid water,
 ! given volFracLiq and soil hydraulic parameters k_sat, theta_sat, theta_res, and m
 ! ******************************************************************************************************************************
 function dHydCond_dLiq(volFracLiq,k_sat,theta_res,theta_sat,m,lTangent) bind(C, name="dHydCond_dLiq")
 import c_double, c_bool
 implicit none
 real(c_double),value :: volFracLiq     ! volumetric fraction of liquid water (-)
 real(c_double),value :: k_sat          ! saturated hydraulic conductivity (m s-1)
 real(c_double),value :: theta_res      ! soil residual volumetric water content (-)
 real(c_double),value :: theta_sat      ! soil porosity (-)
 real(c_double),value :: m              ! vGn "m" parameter (-)
 logical(c_bool),value :: lTangent       ! method used to compute derivative (.true. = analytical)
 real(c_double)       :: dHydCond_dLiq  ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 end function dHydCond_dLiq


 ! ******************************************************************************************************************************
 ! public function dTheta_dTk: differentiate the freezing curve w.r.t. temperature
 ! ******************************************************************************************************************************
 function dTheta_dTk(Tk,theta_res,theta_sat,alpha,n,m) bind(C, name="dTheta_dTk")
 import c_double
 implicit none
 real(c_double),value :: Tk         ! temperature (K)
 real(c_double),value :: theta_res  ! residual liquid water content (-)
 real(c_double),value :: theta_sat  ! porosity (-)
 real(c_double),value :: alpha      ! vGn scaling parameter (m-1)
 real(c_double),value :: n          ! vGn "n" parameter (-)
 real(c_double),value :: m          ! vGn "m" parameter (-)
 real(c_double)       :: dTheta_dTk ! derivative of the freezing curve w.r.t. temperature (K-1)
 end function dTheta_dTk
end interface

contains

 ! ******************************************************************************************************************************
 ! public subroutine iceImpede: compute the ice impedence factor
 ! ******************************************************************************************************************************
 subroutine iceImpede(volFracIce,f_impede, &            ! input
                      iceImpedeFactor,dIceImpede_dLiq)  ! output
 ! computes the ice impedence factor (separate function, as used multiple times)
 implicit none
 ! input variables
 real(rkind),intent(in)     :: volFracIce        ! volumetric fraction of ice (-)
 real(rkind),intent(in)     :: f_impede          ! ice impedence parameter (-)
 ! output variables
 real(rkind)                :: iceImpedeFactor   ! ice impedence factor (-)
 real(rkind)                :: dIceImpede_dLiq   ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
 ! compute ice impedance factor as a function of volumetric ice content
 iceImpedeFactor = 10._rkind**(-f_impede*volFracIce)
 dIceImpede_dLiq = 0._rkind

 end subroutine iceImpede


 ! ******************************************************************************************************************************
 ! public subroutine dIceImpede_dTemp: compute the derivative in the ice impedence factor w.r.t. temperature
 ! ******************************************************************************************************************************
 subroutine dIceImpede_dTemp(volFracIce,dTheta_dT,f_impede,dIceImpede_dT)
 ! computes the derivative in the ice impedance factor w.r.t. temperature
 implicit none
 ! input variables
 real(rkind),intent(in)     :: volFracIce        ! volumetric fraction of ice (-)
 real(rkind),intent(in)     :: dTheta_dT         ! derivative in volumetric liquid water content w.r.t temperature (K-1)
 real(rkind),intent(in)     :: f_impede          ! ice impedence parameter (-)
 ! output variables
 real(rkind)                :: dIceImpede_dT     ! derivative in the ice impedance factor w.r.t. temperature (K-1)
 ! --
 dIceImpede_dT = log(10._rkind)*f_impede*dTheta_dT*10._rkind**(-f_impede*volFracIce)
 end subroutine dIceImpede_dTemp


 ! ******************************************************************************************************************************
 ! public subroutine: compute the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
 ! ******************************************************************************************************************************
 subroutine liquidHead(&
                       ! input
                       matricHeadTotal                          ,& ! intent(in)    : total water matric potential (m)
                       volFracLiq                               ,& ! intent(in)    : volumetric fraction of liquid water (-)
                       volFracIce                               ,& ! intent(in)    : volumetric fraction of ice (-)
                       vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m,& ! intent(in)    : soil parameters
                       dVolTot_dPsi0                            ,& ! intent(in)    : derivative in the soil water characteristic (m-1)
                       dTheta_dT                                ,& ! intent(in)    : derivative in volumetric total water w.r.t. temperature (K-1)
                       ! output
                       matricHeadLiq                            ,& ! intent(out)   : liquid water matric potential (m)
                       dPsiLiq_dPsi0                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                       dPsiLiq_dTemp                            ,& ! intent(out)   : derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                       err,message)                                ! intent(out)   : error control
 ! computes the liquid water matric potential (and the derivatives w.r.t. total matric potential and temperature)
 implicit none
 ! input
 real(rkind),intent(in)            :: matricHeadTotal                           ! total water matric potential (m)
 real(rkind),intent(in)            :: volFracLiq                                ! volumetric fraction of liquid water (-)
 real(rkind),intent(in)            :: volFracIce                                ! volumetric fraction of ice (-)
 real(rkind),intent(in)            :: vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m ! soil parameters
 real(rkind),intent(in)  ,optional :: dVolTot_dPsi0                             ! derivative in the soil water characteristic (m-1)
 real(rkind),intent(in)  ,optional :: dTheta_dT                                 ! derivative in volumetric total water w.r.t. temperature (K-1)
 ! output
 real(rkind),intent(out)           :: matricHeadLiq                             ! liquid water matric potential (m)
 real(rkind),intent(out) ,optional :: dPsiLiq_dPsi0                             ! derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
 real(rkind),intent(out) ,optional :: dPsiLiq_dTemp                             ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 ! output: error control
 integer(i4b),intent(out)       :: err                                       ! error code
 character(*),intent(out)       :: message                                   ! error message
 ! local
 real(rkind)                       :: xNum,xDen                                 ! temporary variables (numeratir, denominator)
 real(rkind)                       :: effSat                                    ! effective saturation (-)
 real(rkind)                       :: dPsiLiq_dEffSat                           ! derivative in liquid water matric potential w.r.t. effective saturation (m)
 real(rkind)                       :: dEffSat_dTemp                             ! derivative in effective saturation w.r.t. temperature (K-1)
 ! ------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='liquidHead/'

 ! ** partially frozen soil
 if(volFracIce > epsilon(1._rkind) .and. matricHeadTotal < 0._rkind)then  ! check that ice exists and that the soil is unsaturated

  ! -----
  ! - compute liquid water matric potential...
  ! ------------------------------------------

  ! - compute effective saturation
  ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
  ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
  xNum   = volFracLiq - theta_res
  xDen   = theta_sat - volFracIce - theta_res
  effSat = xNum/xDen          ! effective saturation

  ! - matric head associated with liquid water
  matricHeadLiq = matricHead(effSat,vGn_alpha,0._rkind,1._rkind,vGn_n,vGn_m)  ! argument is effective saturation, so theta_res=0 and theta_sat=1

  ! compute derivative in liquid water matric potential w.r.t. effective saturation (m)
  if(present(dPsiLiq_dPsi0).or.present(dPsiLiq_dTemp))then
   dPsiLiq_dEffSat = dPsi_dTheta(effSat,vGn_alpha,0._rkind,1._rkind,vGn_n,vGn_m)
  endif

  ! -----
  ! - compute derivative in the liquid water matric potential w.r.t. the total water matric potential...
  ! ----------------------------------------------------------------------------------------------------

  ! check if the derivative is desired
  if(present(dPsiLiq_dTemp))then

   ! (check required input derivative is present)
   if(.not.present(dVolTot_dPsi0))then
    message=trim(message)//'dVolTot_dPsi0 argument is missing'
    err=20; return
   endif

   ! (compute derivative in the liquid water matric potential w.r.t. the total water matric potential)
   dPsiLiq_dPsi0 = dVolTot_dPsi0*dPsiLiq_dEffSat*xNum/(xDen*xDen)

  endif  ! if dPsiLiq_dTemp is desired

  ! -----
  ! - compute the derivative in the liquid water matric potential w.r.t. temperature...
  ! -----------------------------------------------------------------------------------

  ! check if the derivative is desired
  if(present(dPsiLiq_dTemp))then

   ! (check required input derivative is present)
   if(.not.present(dTheta_dT))then
    message=trim(message)//'dTheta_dT argument is missing'
    err=20; return
   endif

   ! (compute the derivative in the liquid water matric potential w.r.t. temperature)
   dEffSat_dTemp = -dTheta_dT*xNum/(xDen*xDen) + dTheta_dT/xDen
   dPsiLiq_dTemp = dPsiLiq_dEffSat*dEffSat_dTemp

  endif  ! if dPsiLiq_dTemp is desired

 ! ** unfrozen soil
 else   ! (no ice)
  matricHeadLiq = matricHeadTotal
  if(present(dPsiLiq_dTemp)) dPsiLiq_dPsi0 = 1._rkind  ! derivative=1 because values are identical
  if(present(dPsiLiq_dTemp)) dPsiLiq_dTemp = 0._rkind  ! derivative=0 because no impact of temperature for unfrozen conditions
 end if  ! (if ice exists)

 end subroutine liquidHead

 ! ******************************************************************************************************************************
 ! public function hydCondMP_liq: compute the hydraulic conductivity of macropores as a function of liquid water content (m s-1)
 ! ******************************************************************************************************************************
 function hydCondMP_liq(volFracLiq,theta_sat,theta_mp,mpExp,satHydCond_ma,satHydCond_mi)
 ! computes hydraulic conductivity given volFracLiq and soil hydraulic parameters
 !  theta_sat, theta_mp, mpExp, satHydCond_ma, and satHydCond_mi
 implicit none
 ! dummies
 real(rkind),intent(in) :: volFracLiq    ! volumetric liquid water content (-)
 real(rkind),intent(in) :: theta_sat     ! soil porosity (-)
 real(rkind),intent(in) :: theta_mp      ! minimum volumetric liquid water content for macropore flow (-)
 real(rkind),intent(in) :: mpExp         ! empirical exponent in macropore flow equation (-)
 real(rkind),intent(in) :: satHydCond_ma ! saturated hydraulic conductivity for macropores (m s-1)
 real(rkind),intent(in) :: satHydCond_mi ! saturated hydraulic conductivity for micropores (m s-1)
 real(rkind)            :: hydCondMP_liq ! hydraulic conductivity (m s-1)
 ! locals
 real(rkind)            :: theta_e     ! effective soil moisture
 if(volFracLiq > theta_mp)then
  theta_e       = (volFracLiq - theta_mp) / (theta_sat - theta_mp)
  hydCondMP_liq = (satHydCond_ma - satHydCond_mi) * (theta_e**mpExp)
 else
  hydCondMP_liq = 0._rkind
 end if
 !write(*,'(a,4(f9.3,1x),2(e20.10))') 'in soil_utils: theta_mp, theta_sat, volFracLiq, hydCondMP_liq, satHydCond_ma, satHydCond_mi = ', &
 !                                                    theta_mp, theta_sat, volFracLiq, hydCondMP_liq, satHydCond_ma, satHydCond_mi
 end function hydCondMP_liq


 ! ******************************************************************************************************************************
 ! public function RH_soilair: compute relative humidity of air in soil pore space
 ! ******************************************************************************************************************************
 function RH_soilair(matpot,Tk)
 implicit none
 real(rkind),intent(in) :: matpot        ! soil water suction -- matric potential (m)
 real(rkind),intent(in) :: Tk            ! temperature (K)
 real(rkind)            :: RH_soilair    ! relative humidity of air in soil pore space
 ! compute relative humidity (UNITS NOTE: Pa = kg m-1 s-2, so R_wv units = m2 s-2 K-1)
 RH_soilair = exp( (gravity*matpot) / (R_wv*Tk) )
 end function RH_soilair


 ! ******************************************************************************************************************************
 ! public function crit_soilT: compute the critical temperature above which all water is unfrozen
 ! ******************************************************************************************************************************
 function crit_soilT(psi)
 implicit none
 real(rkind),intent(in) :: psi           ! matric head (m)
 real(rkind)            :: crit_soilT    ! critical soil temperature (K)
 crit_soilT = Tfreeze + min(psi,0._rkind)*gravity*Tfreeze/LH_fus
 end function crit_soilT


 ! ******************************************************************************************************************************
 ! public function gammp: compute cumulative probability using the Gamma distribution
 ! ******************************************************************************************************************************
 FUNCTION gammp(a,x)
 IMPLICIT NONE
 real(rkind), INTENT(IN) :: a,x
 real(rkind) :: gammp
 if (x<a+1.0_rkind) then
  gammp=gser(a,x)
 else
  gammp=1.0_rkind-gcf(a,x)
 end if
 END FUNCTION gammp


 ! ******************************************************************************************************************************
 ! private function gcf: continued fraction development of the incomplete Gamma function
 ! ******************************************************************************************************************************
 FUNCTION gcf(a,x,gln)
 IMPLICIT NONE
 real(rkind), INTENT(IN) :: a,x
 real(rkind), OPTIONAL, INTENT(OUT) :: gln
 real(rkind) :: gcf
 INTEGER(I4B), PARAMETER :: ITMAX=100
 real(rkind), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
 INTEGER(I4B) :: i
 real(rkind) :: an,b,c,d,del,h
 if (x == 0.0) then
  gcf=1.0
  RETURN
 end if
 b=x+1.0_rkind-a
 c=1.0_rkind/FPMIN
 d=1.0_rkind/b
 h=d
 do i=1,ITMAX
  an=-i*(i-a)
  b=b+2.0_rkind
  d=an*d+b
  if (abs(d) < FPMIN) d=FPMIN
  c=b+an/c
  if (abs(c) < FPMIN) c=FPMIN
  d=1.0_rkind/d
  del=d*c
  h=h*del
  if (abs(del-1.0_rkind) <= EPS) exit
 end do
 if (i > ITMAX) stop 'a too large, ITMAX too small in gcf'
 if (present(gln)) then
  gln=gammln(a)
  gcf=exp(-x+a*log(x)-gln)*h
 else
  gcf=exp(-x+a*log(x)-gammln(a))*h
 end if
 END FUNCTION gcf


 ! ******************************************************************************************************************************
 ! private function gser: series development of the incomplete Gamma function
 ! ******************************************************************************************************************************
 FUNCTION gser(a,x,gln)
 IMPLICIT NONE
 real(rkind), INTENT(IN) :: a,x
 real(rkind), OPTIONAL, INTENT(OUT) :: gln
 real(rkind) :: gser
 INTEGER(I4B), PARAMETER :: ITMAX=100
 real(rkind), PARAMETER :: EPS=epsilon(x)
 INTEGER(I4B) :: n
 real(rkind) :: ap,del,summ
 if (x == 0.0) then
  gser=0.0
  RETURN
 end if
 ap=a
 summ=1.0_rkind/a
 del=summ
 do n=1,ITMAX
  ap=ap+1.0_rkind
  del=del*x/ap
  summ=summ+del
  if (abs(del) < abs(summ)*EPS) exit
 end do
 if (n > ITMAX) stop 'a too large, ITMAX too small in gser'
 if (present(gln)) then
  gln=gammln(a)
  gser=summ*exp(-x+a*log(x)-gln)
 else
  gser=summ*exp(-x+a*log(x)-gammln(a))
 end if
 END FUNCTION gser


 ! ******************************************************************************************************************************
 ! private function gammln: gamma function
 ! ******************************************************************************************************************************
 FUNCTION gammln(xx)
 USE nr_utility_module,only:arth  ! use to build vectors with regular increments
 IMPLICIT NONE
 real(rkind), INTENT(IN) :: xx
 real(rkind) :: gammln
 real(rkind) :: tmp,x
 real(rkind) :: stp = 2.5066282746310005_rkind
 real(rkind), DIMENSION(6) :: coef = (/76.18009172947146_rkind,&
  -86.50532032941677_rkind,24.01409824083091_rkind,&
  -1.231739572450155_rkind,0.1208650973866179e-2_rkind,&
  -0.5395239384953e-5_rkind/)
 if(xx <= 0._rkind) stop 'xx > 0 in gammln'
 x=xx
 tmp=x+5.5_rkind
 tmp=(x+0.5_rkind)*log(tmp)-tmp
 gammln=tmp+log(stp*(1.000000000190015_rkind+&
  sum(coef(:)/arth(x+1.0_rkind,1.0_rkind,size(coef))))/x)
 END FUNCTION gammln


end module soil_utils_module
