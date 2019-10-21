# ===================================================================
# The contents of this file are dedicated to the public domain.  To
# the extent that dedication to the public domain is not available,
# everyone is granted a worldwide, perpetual, royalty-free,
# non-exclusive license to exercise all rights associated with the
# contents of this file for any purpose whatsoever.
# No rights are reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ===================================================================

'''
Created on Fri Oct 19 22:09:37 2019

@author Jérôme Chanut
@author James Harle

$Last commit on:$
'''

import numpy as np


# Function to define an uniform vertical grid:
def set_uniform_refvgrid(dz, jpk):
    depw = np.arange(0, dz * jpk, dz)
    dept = np.arange(0.5 * dz, dz * jpk, dz)
    e3w = np.ones(jpk)*dz
    e3t = np.ones(jpk)*dz
    return depw, dept, e3t, e3w
# Function to define an uniform vertical grid:
def set_zgr_z(dz, jpk):

   rn_hmin     =   -10.    !  min depth of the ocean (>0) or min number of ocean level (<0)
   rn_isfhmin  =    1.00   !  treshold (m) to discriminate grounding ice to floating ice
   rn_e3zps_min=   25.     !  partial step thickness is set larger than the minimum of
   rn_e3zps_rat=    0.2    !  rn_e3zps_min and rn_e3zps_rat*e3t, with 0<rn_e3zps_rat<1
   rn_rdt      =   300.    !  time step for the dynamics (and tracer if nacc=0)   ==> 5760
   jphgr_msh   =       0               !  type of horizontal mesh
   ppglam0     =  999999.0             !  longitude of first raw and column T-point (jphgr_msh = 1)
   ppgphi0     =  999999.0             ! latitude  of first raw and column T-point (jphgr_msh = 1)
   ppe1_deg    =  999999.0             !  zonal      grid-spacing (degrees)
   ppe2_deg    =  999999.0             !  meridional grid-spacing (degrees)
   ppe1_m      =  999999.0             !  zonal      grid-spacing (degrees)
   ppe2_m      =  999999.0             !  meridional grid-spacing (degrees)
   ppsur       =    -3958.951371276829 !  ORCA r4, r2 and r05 coefficients
   ppa0        =     103.9530096000000 ! (default coefficients)
   ppa1        =     2.415951269000000 !
   ppkth       =     15.35101370000000 !
   ppacr       =        7.0            !
   ppdzmin     =  999999.              !  Minimum vertical spacing
   pphmax      =  999999.              !  Maximum depth
   ldbletanh   =    .TRUE.             !  Use/do not use double tanf function for vertical coordinates
   ppa2        =      100.760928500000 !  Double tanh function parameters
   ppkth2      =       48.029893720000 !
   ppacr2      =       13.000000000000 !



    IF( nn_timing == 1 )  CALL timing_start('zgr_z')
      !
      ! Set variables from parameters
      ! ------------------------------
       zkth = ppkth       ;   zacr = ppacr
       zdzmin = ppdzmin   ;   zhmax = pphmax
       zkth2 = ppkth2     ;   zacr2 = ppacr2   ! optional (ldbletanh=T) double tanh parameters

      ! If ppa1 and ppa0 and ppsur are et to pp_to_be_computed
      !  za0, za1, zsur are computed from ppdzmin , pphmax, ppkth, ppacr
      IF(   ppa1  == pp_to_be_computed  .AND.  &
         &  ppa0  == pp_to_be_computed  .AND.  &
         &  ppsur == pp_to_be_computed           ) THEN
         !
         za1  = (  ppdzmin - pphmax / FLOAT(jpkm1)  )                                                      &
            & / ( TANH((1-ppkth)/ppacr) - ppacr/FLOAT(jpk-1) * (  LOG( COSH( (jpk - ppkth) / ppacr) )      &
            &                                                   - LOG( COSH( ( 1  - ppkth) / ppacr) )  )  )
         za0  = ppdzmin - za1 *              TANH( (1-ppkth) / ppacr )
         zsur =   - za0 - za1 * ppacr * LOG( COSH( (1-ppkth) / ppacr )  )
      ELSE
         za1 = ppa1 ;       za0 = ppa0 ;          zsur = ppsur
         za2 = ppa2                            ! optional (ldbletanh=T) double tanh parameter
      ENDIF

      IF(lwp) THEN                         ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates'
         WRITE(numout,*) '    ~~~~~~~'
         IF(  ppkth == 0._wp ) THEN
              WRITE(numout,*) '            Uniform grid with ',jpk-1,' layers'
              WRITE(numout,*) '            Total depth    :', zhmax
              WRITE(numout,*) '            Layer thickness:', zhmax/(jpk-1)
         ELSE
            IF( ppa1 == 0._wp .AND. ppa0 == 0._wp .AND. ppsur == 0._wp ) THEN
               WRITE(numout,*) '         zsur, za0, za1 computed from '
               WRITE(numout,*) '                 zdzmin = ', zdzmin
               WRITE(numout,*) '                 zhmax  = ', zhmax
            ENDIF
            WRITE(numout,*) '           Value of coefficients for vertical mesh:'
            WRITE(numout,*) '                 zsur = ', zsur
            WRITE(numout,*) '                 za0  = ', za0
            WRITE(numout,*) '                 za1  = ', za1
            WRITE(numout,*) '                 zkth = ', zkth
            WRITE(numout,*) '                 zacr = ', zacr
            IF( ldbletanh ) THEN
               WRITE(numout,*) ' (Double tanh    za2  = ', za2
               WRITE(numout,*) '  parameters)    zkth2= ', zkth2
               WRITE(numout,*) '                 zacr2= ', zacr2
            ENDIF
         ENDIF
      ENDIF

      ! Reference z-coordinate (depth - scale factor at T- and W-points)
      ! ======================
      IF( ppkth == 0._wp ) THEN            !  uniform vertical grid 



         za1 = zhmax / FLOAT(jpk-1)

         DO jk = 1, jpk
            zw = FLOAT( jk )
            zt = FLOAT( jk ) + 0.5_wp
            gdepw_1d(jk) = ( zw - 1 ) * za1
            gdept_1d(jk) = ( zt - 1 ) * za1
            e3w_1d  (jk) =  za1
            e3t_1d  (jk) =  za1
         END DO
      ELSE                                ! Madec & Imbard 1996 function
         IF( .NOT. ldbletanh ) THEN
            DO jk = 1, jpk
               zw = REAL( jk , wp )
               zt = REAL( jk , wp ) + 0.5_wp
               gdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr * LOG ( COSH( (zw-zkth) / zacr ) )  )
               gdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr * LOG ( COSH( (zt-zkth) / zacr ) )  )
               e3w_1d  (jk) =          za0      + za1        * TANH(       (zw-zkth) / zacr   )
               e3t_1d  (jk) =          za0      + za1        * TANH(       (zt-zkth) / zacr   )
            END DO
         ELSE
            DO jk = 1, jpk
               zw = FLOAT( jk )
               zt = FLOAT( jk ) + 0.5_wp
               ! Double tanh function
               gdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr * LOG ( COSH( (zw-zkth ) / zacr  ) )    &
                  &                             + za2 * zacr2* LOG ( COSH( (zw-zkth2) / zacr2 ) )  )
               gdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr * LOG ( COSH( (zt-zkth ) / zacr  ) )    &
                  &                             + za2 * zacr2* LOG ( COSH( (zt-zkth2) / zacr2 ) )  )
               e3w_1d  (jk) =          za0      + za1        * TANH(       (zw-zkth ) / zacr  )      &
                  &                             + za2        * TANH(       (zw-zkth2) / zacr2 )
               e3t_1d  (jk) =          za0      + za1        * TANH(       (zt-zkth ) / zacr  )      &
                  &                             + za2        * TANH(       (zt-zkth2) / zacr2 )
            END DO
         ENDIF
         gdepw_1d(1) = 0._wp                    ! force first w-level to be exactly at zero
      ENDIF

      IF ( ln_isfcav .OR. ln_e3_dep ) THEN      ! e3. = dk[gdep]   
         !
!==>>>   need to be like this to compute the pressure gradient with ISF. 
!        If not, level beneath the ISF are not aligned (sum(e3t) /= depth)
!        define e3t_0 and e3w_0 as the differences between gdept and gdepw respectively
!
         DO jk = 1, jpkm1
            e3t_1d(jk) = gdepw_1d(jk+1)-gdepw_1d(jk)
         END DO
         e3t_1d(jpk) = e3t_1d(jpk-1)   ! we don't care because this level is masked in NEMO

         DO jk = 2, jpk
            e3w_1d(jk) = gdept_1d(jk) - gdept_1d(jk-1)
         END DO
         e3w_1d(1  ) = 2._wp * (gdept_1d(1) - gdepw_1d(1))
      END IF


    
    depw = np.arange(0, dz * jpk, dz)
    dept = np.arange(0.5 * dz, dz * jpk, dz)
    e3w = np.ones(jpk)*dz
    e3t = np.ones(jpk)*dz
    return depw, dept, e3t, e3w


# Function to define vertical grid parameters in the full cells case:
def set_zcovgrid(bat, depw_1d, dept_1d, e3t_1d, e3w_1d):
    (jpi, jpj) = np.shape(bat)
    jpk = np.size(depw_1d)
    # set bottom level
    kbot = np.ones((jpi, jpj), dtype=int)*jpk - 1
    for k in np.arange(jpk - 2, -1, -1):
        kbot = np.where((bat < dept_1d[k]), k, kbot)

    # set scale factors and depths at T-points:
    e3t = np.zeros((jpi, jpj, jpk))
    e3w = np.zeros((jpi, jpj, jpk))
    depw = np.zeros((jpi, jpj, jpk))
    dept = np.zeros((jpi, jpj, jpk))
    bato = np.zeros((jpi, jpj))
    for k in np.arange(0, jpk, 1):
        e3t[:, :, k] = e3t_1d[k]
        e3w[:, :, k] = e3w_1d[k]
        depw[:, :, k] = depw_1d[k]
        dept[:, :, k] = dept_1d[k]

    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            for k in np.arange(0, kbot[i, j], 1):
                bato[i, j] = bato[i, j] + e3t[i, j, k]

    return kbot, bato, e3t, e3w, depw, dept


# Function to define vertical grid parameters in partial cell case:
def set_pstepvgrid(bat, depw_1d, dept_1d, e3t_1d, e3w_1d):
    (jpi, jpj) = np.shape(bat)
    jpk = np.size(depw_1d)
    # set bottom level
    kbot = np.ones((jpi, jpj), dtype=int)*jpk - 1
    for k in np.arange(jpk - 2, -1, -1):
        zmin = 0.1 * e3t_1d[k]
        kbot = np.where((bat < (depw_1d[k] + zmin)), k, kbot)

    # set scale factors and depths at T-points:
    e3t = np.zeros((jpi, jpj, jpk))
    e3w = np.zeros((jpi, jpj, jpk))
    depw = np.zeros((jpi, jpj, jpk))
    dept = np.zeros((jpi, jpj, jpk))
    bato = np.zeros((jpi, jpj))
    for k in np.arange(0, jpk, 1):
        e3t[:, :, k] = e3t_1d[k]
        e3w[:, :, k] = e3w_1d[k]
        depw[:, :, k] = depw_1d[k]
        dept[:, :, k] = dept_1d[k]

    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            k = np.int(kbot[i, j]) - 1
            if k >= 0:
                depw[i, j, k + 1] = min(bat[i, j], depw_1d[k + 1])
                e3t[i, j, k] = depw[i, j, k + 1] - depw[i, j, k]
                dept[i, j, k] = depw[i, j, k] + 0.5 * e3t[i, j, k]
                e3w[i, j, k] = dept[i, j, k] - dept[i, j, k - 1]

    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            for k in np.arange(0, kbot[i, j], 1):
                bato[i, j] = bato[i, j] + e3t[i, j, k]

    return kbot, bato, e3t, e3w, depw, dept


# Function to define vertical grid parameters in the s-coordinate case:
def set_scovgrid(bat, depw_1d, dept_1d, e3t_1d, e3w_1d):
    (jpi, jpj) = np.shape(bat)
    jpk = np.size(depw_1d)
    jpkm1 = jpk - 1
    # set bottom level
    kbot = np.ones((jpi, jpj)) * jpkm1

    # set scale factors and depths at T-points:
    e3t = np.zeros((jpi, jpj, jpk))
    e3w = np.zeros((jpi, jpj, jpk))
    depw = np.zeros((jpi, jpj, jpk))
    dept = np.zeros((jpi, jpj, jpk))
    for k in np.arange(0, jpk, 1):
        e3t[:, :, k] = e3t_1d[k]
        e3w[:, :, k] = e3w_1d[k]
        depw[:, :, k] = depw_1d[k]
        dept[:, :, k] = dept_1d[k]

    # Uniform sigma distribution:
    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            if bat[i, j] > 0:
                e3t[i, j, :] = bat[i, j] / np.float(jpkm1)
                e3w[i, j, :] = bat[i, j] / np.float(jpkm1)
                e3w[i, j, 0] = 0.5 * bat[i, j] / np.float(jpkm1)
            dept[i, j, :] = np.cumsum(e3w[i, j, :])
            depw[i, j, :] = np.cumsum(e3t[i, j, :]) - e3t[i, j, 0]

    return kbot, bat, e3t, e3w, depw, dept


# Function to define vertical grid parameters in the s-coordinate on top of z case:
def set_scotopofzvgrid(bat, depmax, depw_1d, dept_1d, e3t_1d, e3w_1d):
    (jpi, jpj) = np.shape(bat)
    jpk = np.size(depw_1d)
    jpkm1 = jpk - 1

    # Set number of s levels:
    jpks = 0
    for k in range(0, jpk-1):
        if depw_1d[k+1] <= depmax:
            jpks = k

    # set bottom level in the z case first
    kbot = np.ones((jpi, jpj), dtype=int) * jpkm1
    for k in np.arange(jpk - 2, -1, -1):
        kbot = np.where((bat < dept_1d[k]), k, kbot)

    # set s coordinates on top:
    kbot = np.where(kbot < jpks, jpks, kbot)

    # set scale factors and depths at T-points:
    e3t = np.zeros((jpi, jpj, jpk))
    e3w = np.zeros((jpi, jpj, jpk))
    depw = np.zeros((jpi, jpj, jpk))
    dept = np.zeros((jpi, jpj, jpk))
    for k in np.arange(0, jpk, 1):
        e3t[:, :, k] = e3t_1d[k]
        e3w[:, :, k] = e3w_1d[k]
        depw[:, :, k] = depw_1d[k]
        dept[:, :, k] = dept_1d[k]

    # Uniform sigma distribution for the top jpks levels:
    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            if bat[i, j] > 0:
                zbat = min(bat[i, j], depw_1d[jpks])
                for k in range(0, jpks+1):
                    e3t[i, j, k] = zbat * e3t_1d[k] / depw_1d[jpks]
                    e3w[i, j, k] = zbat * e3w_1d[k] / depw_1d[jpks]
                e3w[i, j, 0] = 0.5 * e3t[i, j, 0]

    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            if bat[i, j] > 0.:
                dept[i, j, :] = np.cumsum(e3w[i, j, :])
                depw[i, j, :] = np.cumsum(e3t[i, j, :]) - e3t[i, j, 0]

    # Update bathymetry:
    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            bat[i, j] = depw[i,j,kbot[i,j]]

    return kbot, bat, e3t, e3w, depw, dept


def set_uvfvgrid(ln_zps, ln_sco, e3t_1d, e3w_1d, e3t, e3w):
    #
    (jpi, jpj, jpk) = np.shape(e3t)
    e3u = np.zeros((jpi, jpj, jpk))
    e3v = np.zeros((jpi, jpj, jpk))
    e3uw = np.zeros((jpi, jpj, jpk))
    e3vw = np.zeros((jpi, jpj, jpk))
    e3f = np.zeros((jpi, jpj, jpk))

    for k in np.arange(0, jpk, 1):
        e3u[:, :, k] = e3t_1d[k]
        e3v[:, :, k] = e3t_1d[k]
        e3f[:, :, k] = e3t_1d[k]
        e3uw[:, :, k] = e3w_1d[k]
        e3vw[:, :, k] = e3w_1d[k]

    if ln_zps == 1:
        for i in np.arange(0, jpi - 1, 1):
            for j in np.arange(0, jpj - 1, 1):
                for k in np.arange(0, jpk - 1, 1):
                    e3u[i, j, k] = min(e3t[i, j, k], e3t[i + 1, j, k])
                    e3v[i, j, k] = min(e3t[i, j, k], e3t[i, j + 1, k])
                    e3uw[i, j, k] = min(e3w[i, j, k], e3w[i + 1, j, k])
                    e3vw[i, j, k] = min(e3w[i, j, k], e3w[i, j + 1, k])

        for i in np.arange(0, jpi - 1, 1):
            for j in np.arange(0, jpj - 1, 1):
                for k in np.arange(0, jpk - 1, 1):
                    e3f[i, j, k] = min(e3v[i, j, k], e3v[i + 1, j, k])
    elif ln_sco == 1:
        for i in np.arange(0, jpi - 1, 1):
            for j in np.arange(0, jpj - 1, 1):
                for k in np.arange(0, jpk - 1, 1):
                    e3u[i, j, k] = 0.5 * (e3t[i, j, k] + e3t[i + 1, j, k])
                    e3v[i, j, k] = 0.5 * (e3t[i, j, k] + e3t[i, j + 1, k])
                    e3uw[i, j, k] = 0.5 * (e3w[i, j, k] + e3w[i + 1, j, k])
                    e3vw[i, j, k] = 0.5 * (e3w[i, j, k] + e3w[i, j + 1, k])

        for i in np.arange(0, jpi - 1, 1):
            for j in np.arange(0, jpj - 1, 1):
                for k in np.arange(0, jpk - 1, 1):
                    e3f[i, j, k] = 0.5 * (e3v[i, j, k] + e3v[i + 1, j, k])

    return e3u, e3v, e3f, e3uw, e3vw


def set_scovgrid_step(bat, depw_1d, dept_1d, e3t_1d, e3w_1d):
    (jpi, jpj) = np.shape(bat)
    jpk = np.size(depw_1d)
    jpkm1 = jpk - 1
    # set bottom level
    kbot = np.ones((jpi, jpj)) * jpkm1
    # set scale factors and depths at T-points:
    e3t = np.zeros((jpi, jpj, jpk))
    e3w = np.zeros((jpi, jpj, jpk))
    depw = np.zeros((jpi, jpj, jpk))
    dept = np.zeros((jpi, jpj, jpk))
    for k in np.arange(0, jpk, 1):
        e3t[:, :, k] = e3t_1d[k]
        e3w[:, :, k] = e3w_1d[k]
        depw[:, :, k] = depw_1d[k]
        dept[:, :, k] = dept_1d[k]

    # Uniform sigma distribution:
    for i in np.arange(0, jpi, 1):
        for j in np.arange(0, jpj, 1):
            if bat[i, j] > 0:
                e3t[i, j, :] = bat[i, j] / np.float64(jpkm1)
                e3w[i, j, :] = bat[i, j] / np.float64(jpkm1)
                e3w[i, j, 0] = 0.5 * bat[i, j] / np.float64(jpkm1)
                #
                dept[i, j, :] = np.cumsum(e3w[i, j, :])
                depw[i, j, :] = np.cumsum(e3t[i, j, :]) - e3t[i, j, 0]

    e3u = np.zeros((jpi, jpj, jpk))
    e3v = np.zeros((jpi, jpj, jpk))
    e3uw = np.zeros((jpi, jpj, jpk))
    e3vw = np.zeros((jpi, jpj, jpk))
    e3f = np.zeros((jpi, jpj, jpk))

    for k in np.arange(0, jpk, 1):
        e3u[:, :, k] = e3t_1d[k]
        e3v[:, :, k] = e3t_1d[k]
        e3f[:, :, k] = e3t_1d[k]
        e3uw[:, :, k] = e3w_1d[k]
        e3vw[:, :, k] = e3w_1d[k]

    batus = np.zeros((jpi, jpj))
    batvs = np.zeros((jpi, jpj))
    batu = np.zeros((jpi, jpj))
    batv = np.zeros((jpi, jpj))
    for i in np.arange(0, jpi - 1, 1):
        for j in np.arange(0, jpj - 1, 1):
            batus[i, j] = min(bat[i, j], bat[i + 1, j])
            batu[i, j] = 0.5 * (bat[i, j] + bat[i + 1, j])
            batvs[i, j] = min(bat[i, j], bat[i, j + 1])
            batv[i, j] = 0.5 * (bat[i, j] + bat[i, j + 1])

    for i in np.arange(0, jpi - 1, 1):
        for j in np.arange(0, jpj - 1, 1):
            zbot = batus[i, j]
            for k in np.arange(jpk - 2, -1, -1):
                zsup = 0.5 * (depw[i, j, k] + depw[i + 1, j, k])
                if batus[i, j] < batu[i, j]:    # step to the right
                    zsup = min(zsup, zbot - 0.8 * e3t[i + 1, j, k])
                else:                           # step to the left
                    zsup = min(zsup, zbot - 0.8 * e3t[i, j, k])
                e3u[i, j, k] = zbot - zsup
                zbot = zsup

    return kbot, bat, e3t, e3w, depw, dept, e3u
