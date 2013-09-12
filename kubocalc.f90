program kubocalc
!
! optical conductivity via Kubo-Greenwood using pwscf output files
!
! Lazaro Calderin,  calderin@lcalderin.net
! Jan Beck       ,  kubocalc@janbeck.com

! important references:
! - Calderin et al.; Structural, dynamic, and electronic properties of liquid tin: An ab initio molecular dynamics study
!   For the form of Kubo-Greenwood equation used for sigma/Lxx
!
! - Holst et al.   ; Electronic transport coefﬁcients from ab initio simulations and application to dense liquid hydrogen
!   For the relationship of sigma/Lxx to thermopower and electronic thermal conductivity
!   note that one unit of charge is absorbed into the formula for L12. Moving that back out, one gets the more commonly seen
!   S = 1/eT (L12/L11)   

! Global variables
!

use kinds,      ONLY : DP
use io_files,   ONLY : nd_nmbr, prefix, outdir, tmp_dir, nwordwfc, iunwfc
use io_global,  ONLY : ionode
use ener,       ONLY : ef
use klist,      ONLY : nks, xk, wk, nelec, ngauss
use constants,  ONLY : PI, RYTOEV, tpi, ELECTRON_SI, H_PLANCK_SI, BOHR_RADIUS_SI, ELECTRONMASS_SI, eps6
use cell_base,  ONLY : tpiba, tpiba2, omega
use wvfct,      ONLY : npw, npwx, nbnd, igk, g2kin, wg, et
use gvect,      ONLY : ngm, g, ecutwfc
use lsda_mod,   ONLY : lsda, nspin
USE uspp,       ONLY : okvan
use wavefunctions_module, ONLY : evc
use input_parameters, ONLY : calculation, restart_mode, pseudo_dir, disk_io, wfcdir, wf_collect
Use noncollin_module, ONLY : noncolin, npol
Use environment, ONLY : environment_start

implicit none
 

  !
  ! Local variables
  !

  logical  :: check_wfc=.false.
  real(DP) :: temperature=800.0_DP
  character(len=4) :: property='all', deltarep='l'

  integer  :: ios, ierr ! Input and memory allocation status variables

  integer  :: nw        ! Number of frequency points
  real(DP) :: wmax, dw  ! Maximum frequency and frequency step
  real(DP) :: sg, ng    ! Gaussian parameters
  real(DP) :: sa        ! Lorentzian width  
  integer  :: iw        ! Frequency index
  real(DP), allocatable    :: wgrid(:) ! Frequency grid
  real(DP), allocatable    :: fermi_distribution(:) ! Fermi distribution at current temp
  real(DP), allocatable    :: fermi_distribution_derivative_neg(:) ! derivative of Fermi distribution at current temp
  real(DP)                 :: conductivity_sum        ! for integral of conductivity trace

  integer  :: ik, ig, ib, ib1, ib2, ix1, ix2 ! ik: k-point index, ig: plane wave index, ib: band index
                                             ! ix: cartesian index

  real(DP), allocatable    :: focc(:,:)  ! State occupations (fermi-dirac)

  complex(DP), allocatable    :: dipolem(:,:,:,:), sigma(:,:,:)      ! Dipole matrix elements andconductivity
  complex(DP), allocatable    :: L11(:,:,:),L12(:,:,:),L21(:,:,:), L12_symmetric(:,:,:),L22aa(:,:,:) ,L22ab(:,:,:)  ,L22bb(:,:,:) 
  complex(DP), allocatable    :: thermopower(:,:,:),thermopower_symmetric(:,:,:),thermopower1(:,:,:),thermopower2(:,:,:)  
  complex(DP), allocatable    :: thermal_conductivity(:,:,:)

  complex(DP) :: dipole_aux(3),  thesum,thesum1,thesum2 ! auxiliar variables
  complex(DP) :: caux, cmaux(3,3), cmaux1(3,3),cmaux2(3,3),cmaux22aa(3,3),cmaux22ab(3,3),cmaux22bb(3,3),cmaux_symmetric(3,3)
  real(DP)    :: aux, ee, e_minus_mu1, e_minus_mu2 , e_minus_mu_symmetric     


  NAMELIST / control / calculation, restart_mode, outdir, pseudo_dir, prefix, disk_io, wfcdir, wf_collect
  NAMELIST / inputpp / check_wfc, property, deltarep, temperature

  !
  ! constant definitions
  !
  nw   = 1000          ! Number of points in the frequency grid
  wmax = 6.0_dp        ! Maximum frequency in eV
  dw   = 0.001         ! Frecuency step
  nw   = int(wmax/dw)
  write(*,*) "Number of points in the frequency grid: ", nw

  call environment_start( 'POST-PROC' )
  !
  ! input
  !

  call input_from_file( )

  if ( ionode )  then 
     !
     read (5, control, IOSTAT=ios)
     if (ios/=0) call errore('optcond', 'reading namelist CONTROL', ABS(ios))
     !
     tmp_dir = TRIM(outdir)
     read(5, inputpp, IOSTAT=ios)
     if( ios/=0 ) call errore('optcond', 'reading namelist INPUTPP', ABS(ios))
  endif

  call read_file
  call openfil_pp

  write(*,*) ' '   

  !
  ! Check if this program can be used
  !

  if( okvan )                 call errore( 'optcond','ultrasoft not implemented'    , 1 )
  if( noncolin )              call errore( 'optcond','non-collinear not implemented', 2 )
  if( wf_collect .eqv. .false. ) call errore( 'optcond','wave functions not collected' , 3 ) 

  !
  ! Check the kind of smearing used by pw.x
  !

  select case (ngauss)
  case ( -99 )
  case default
    write(*,*) "Warning: not fermi-dirac smearing, ngauss=", ngauss
  end select

  !
  ! lets make some space ...
  !
  allocate( sigma(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating sigma(w)', ABS(ierr) )
  allocate( L11(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating L11(w)', ABS(ierr) )
  allocate( L12(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating L12(w)', ABS(ierr) )
  allocate( L21(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating L21(w)', ABS(ierr) )
  allocate( L22aa(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating L22aa(w)', ABS(ierr) )
  allocate( L22ab(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating L22ab(w)', ABS(ierr) )
  allocate( L22bb(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating L22bb(w)', ABS(ierr) )
  allocate( L12_symmetric(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating L12_symmetric(w)', ABS(ierr) )
  allocate( thermopower(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating thermopower(w)', ABS(ierr) )
  allocate( thermal_conductivity(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating thermal_conductivity(w)', ABS(ierr) )
  allocate( thermopower_symmetric(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating thermopower_symmetric(w)', ABS(ierr) )
  allocate( thermopower1(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating thermopower1(w)', ABS(ierr) )
  allocate( thermopower2(3,3,nw), STAT=ierr )
  if (ierr/=0) CALL errore('optcond','allocating thermopower2(w)', ABS(ierr) )
  allocate(focc(nbnd,nks), STAT=ierr)                                           ! Ocuppation numbers (separating the fermi-dirac)
  if (ierr/=0) call errore('optcond','allocating focc',ABS(ierr))
  allocate( dipolem(nks,3,nbnd, nbnd), STAT=ierr )                              ! Dipole matrix elements 
      if (ierr/=0) call errore('optcond','allocating dipole', abs(ierr) )


  ! according to pwscf wg=wk*fermidirac
 
  focc=0.0_dp
  do ik = 1, nks
     do ib  = 1,nbnd
        focc(ib,ik)= wg(ib,ik)*(3-nspin)/wk(ik)
     enddo
  enddo

  write(*,*) '--------------------------------------------'
  write(*,*) 'Number of electrons                          :', nelec
  write(*,*) 'Number of bands                              :', nbnd
  write(*,*) 'Number of k-points                           :', nks
  write(*,*) 'Sum of k-weights                             :', sum(wk)         ! Should be equal to 2
  write(*,*) 'Sum of band weights                          :', sum(wg)         ! Should be equal to nelec
  write(*,*) 'Sum of fermi-occup                           :', sum(focc(:,:))  ! Should be equal to the total number of occupied states
  write(*,*) 'Volume (bohr^3)                              :', omega
  write(*,*) 'The Fermi energy is (eV)                     :', ef*rytoev
  write(*,*) 'The temperature for the Fermi distribution is:', temperature

  !
  ! Checking ortho-normalization of the eigenvectors
  !

!$OMP PARALLEL
write(*,*) "Hello. Starting an OMP thread...."
!$OMP END PARALLEL


  select case ( check_wfc )
  case (.true.)

     do ik = 1, nks

        write(*,*) 'Checking orthonormalization of bands for the kpoint:',ik 

        ! Get npw for the current k-vector
        call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)

        ! Read the wave function for the current k-vector
        call davcio (evc , nwordwfc, iunwfc, ik, - 1)

        do ib1 = 1, nbnd
           do ib2 = 1, nbnd
           
              ! Get the g-vectors of the current k-vector
              call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)  

              ! Read the wave function 
              call davcio (evc , nwordwfc, iunwfc, ik, - 1)                          

              caux=0.0_dp

              do ig=1, npw
                 caux= caux + conjg( evc(ig,ib1) ) * evc(ig,ib2)
              enddo
        
              if ( (ib1/=ib2) .and. (abs(real(caux,dp)) > 1.d-3 .or. abs(aimag(caux))> 1.d-3) ) then
                 write(*,*) 'Non orthogonal:', ik, ib1, ib2, caux
              endif 

              if ( (ib1==ib2) .and. abs(real(caux,dp)-1.0_dp)>1.0d-3 ) then
                 write(*,*) 'Non normalized:', ik, ib1, ib2, caux
              endif

           enddo
        enddo   

     enddo

  end select

  !
  ! Dipole matrix elements 
  ! 

  write(*,*) ''
  write(*,*) 'Calculating the matrix elements ..'

  dipolem = ( 0.0_dp, 0.0_dp )

  do ik = 1, nks  
     write(*,'(a,i4)') 'Calculating the matrix elements for k-point number:', ik

     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin) ! Get the G vectors for the current k-vector  
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)                         ! Read the wave function for the given k-vector 
     !
     ! Calculate dipole matrix elements
     !
     !$OMP DO ORDERED
     do ib1=1, nbnd
        do ib2=1, nbnd        
           !$OMP ORDERED   
           dipole_aux = ( 0.0_dp, 0.0_dp )
           do  ig=1,npw
               caux= conjg( evc(ig,ib1) ) * evc(ig,ib2) 
               dipole_aux(:) = dipole_aux(:) + ( xk(:,ik) + g(:,igk(ig)) ) * caux
           enddo
         
           dipolem(ik,:,ib1,ib2) = tpiba * dipole_aux(:)
           !$OMP END ORDERED
        enddo
     enddo
     !$OMP END DO

  enddo

  !
  ! Frequency grid
  !

  aux = et(2,1) - et(1,1)
  do ik=2, nks
        if( (et(2,ik)-et(1,ik)) > aux ) aux = et(2,ik)-et(1,ik) 
  enddo

 
  allocate(wgrid(nw+1),STAT=ierr)
  if (ierr/=0) call errore('optcond','allocating wgrid', abs(ierr))
  allocate(fermi_distribution(nw+1),STAT=ierr)
  if (ierr/=0) call errore('optcond','allocating fermi_distribution', abs(ierr))
  allocate(fermi_distribution_derivative_neg(nw+1),STAT=ierr)
  if (ierr/=0) call errore('optcond','allocating fermi_distribution_derivative_neg', abs(ierr))
  do iw = 1, nw
     wgrid(iw) = iw*dw
     fermi_distribution(iw)=1.0_DP/(EXP(wgrid(iw)/(0.00008617332_DP*temperature))+1.0_DP)
     if (iw .eq. 1) then
         fermi_distribution_derivative_neg(iw) = (0.5_DP-fermi_distribution(iw))/dw
     else
         fermi_distribution_derivative_neg(iw) = (fermi_distribution(iw-1)-fermi_distribution(iw))/dw
     end if
  enddo

  write(*,*) "Frequency step (eV):", dw

  
  write(*,'(a)',advance='no') 'Calculating the optical conductivity ..'

  select case ( deltarep )
  case ( 'l' )

     sa = 0.5d0*dw ! smearing parameter in the Lorenztian (for the Harrison version)
     !sa = 0.2_dp
     write(*,'(a)') ' Using a Lorentzian for the delta function'
     write(*,*)     '    Lorentzian width (broadening parameter) in eV:', sa
     !$OMP DO ORDERED
     do iw = 1, nw
        !$OMP ORDERED 
        cmaux  = (0.0_DP,0.0_DP)
        cmaux1 = (0.0_DP,0.0_DP)
        cmaux2 = (0.0_DP,0.0_DP)
        cmaux_symmetric = (0.0_DP,0.0_DP)
        write(*,'(a,i10,a,i10)') 'Calculating frequency grid point number:', iw, ' out of ', nw
        do ik = 1, nks

           do ib1=1, nbnd
              !if( focc(ib1,ik) < eps6 ) cycle

              do ib2=1, nbnd
                 !if( (1.0_dp-eps6) < focc(ib2,ik) ) cycle

                 ee = ( et(ib2,ik) - et(ib1,ik) )*RYTOEV - wgrid(iw)
                 e_minus_mu1 = (et(ib1,ik) - ef) * RYTOEV
                 e_minus_mu2 = (et(ib2,ik) - ef) * RYTOEV
                 e_minus_mu_symmetric = ( (et(ib1,ik) + et(ib2,ik)) / 2 - ef) * RYTOEV
                 do ix1 = 1, 3
                    do ix2 = 1, 3 
                       thesum = ( wg(ib1,ik) - wg(ib2,ik) )   
                       thesum1 = conjg( dipolem(ik,ix1,ib1,ib2) )   
                       thesum2 = dipolem(ik,ix2,ib1,ib2) * sa /(ee*ee+sa*sa)/ wgrid(iw)/PI
                       cmaux(ix1, ix2)           = cmaux(ix1,ix2)                                       + thesum2 * thesum1 * thesum
                       cmaux1(ix1, ix2)          = cmaux1(ix1,ix2)          + e_minus_mu1               * thesum2 * thesum1 * thesum
                       cmaux2(ix1, ix2)          = cmaux2(ix1,ix2)          + e_minus_mu2               * thesum2 * thesum1 * thesum
                       cmaux_symmetric(ix1, ix2) = cmaux_symmetric(ix1,ix2) + e_minus_mu_symmetric      * thesum2 * thesum1 * thesum
                       cmaux22aa(ix1, ix2)       = cmaux22aa(ix1,ix2)       + e_minus_mu1 * e_minus_mu2 * thesum2 * thesum1 * thesum
                       cmaux22ab(ix1, ix2)       = cmaux22ab(ix1,ix2)       + e_minus_mu1 * e_minus_mu2 * thesum2 * thesum1 * thesum
                       cmaux22bb(ix1, ix2)       = cmaux22bb(ix1,ix2)       + e_minus_mu_symmetric*e_minus_mu_symmetric * thesum2 * thesum1 * thesum
                    enddo
                 enddo

              enddo
           enddo

        enddo
       
       sigma(:,:,iw)           = cmaux(:,:)
       L12(:,:,iw)             =  cmaux1(:,:)
       L21(:,:,iw)             =  cmaux2(:,:) 
       L12_symmetric(:,:,iw)   =  cmaux_symmetric(:,:) 
       L22aa(:,:,iw)           =  cmaux22aa(:,:) 
       L22ab(:,:,iw)           =  cmaux22ab(:,:) 
       L22bb(:,:,iw)           =  cmaux22bb(:,:) 
       !$OMP END ORDERED
     enddo
     !$OMP END DO
  case ( 'g' )

     sg  = 10.0d0 * dw                        ! Gaussian width in eV
     ng  =  1.0d0 / ( sg*sqrt(2.0d0*PI) )     ! Gaussian normalization factor
     sg  =  1.0d0 / ( 2.0d0*sg*sg )           ! using sg as an aux variable for the exponent
                                              ! of the Gaussian

     write(*,'(a)') ' using a Gaussian for the delta function' 
     write(*,*) 'Gaussian width (eV):', sa

     do iw = 1, nw
        cmaux = (0.0d0,0.0d0)
   
        do ik = 1, nks

           do ib1=1, nbnd
              !if( focc(ib1,ik) < eps6 ) cycle

              do ib2=1, nbnd
                 !if( (1.0_dp-eps6) < focc(ib2,ik) ) cycle

                 ee = ( et(ib2,ik) - et(ib1,ik) )*RYTOEV - wgrid(iw)

                 do ix1 = 1, 3
                    do ix2 = 1, 3
                       thesum = ( wg(ib1,ik) - wg(ib2,ik) )
                       thesum = thesum * conjg( dipolem(ik,ix1,ib1,ib2) ) 
                       thesum = thesum * dipolem(ik,ix2,ib1,ib2) 
                       thesum = thesum * ng * exp(-ee*ee*sg ) / wgrid(iw)
                       cmaux(ix1, ix2) = cmaux(ix1,ix2) + thesum
                    enddo
                 enddo

              enddo
           enddo

        enddo

        sigma(:,:,iw) = cmaux(:,:)

     enddo

  case default
     write(*,*) 'Representation of dirac delta function not implemented:', deltarep

  end select  

  !
  ! Setting up the constant 2 Pi electron_charge^2/electron_mass^2 1/volume_unit_cell
  ! in SI units and also the change of units to get the units of 1/(ohm*meter) for the 
  ! conductivity
  !

  aux = PI * ELECTRON_SI**2 * (H_PLANCK_SI/tpi)**2 / & 
        (ELECTRONMASS_SI**2*omega*BOHR_RADIUS_SI**3)  

  aux = aux * H_PLANCK_SI/tpi * 6.24150974d18            ! Correcting for w given as hbar w and in eV
  aux = aux / BOHR_RADIUS_SI**2                          ! Correcting for grad. given in 1/Bohr^2
  aux = aux * 6.24150974d18                              ! Correcting for sg given in eV  
  aux = aux / 1.0E8_DP                                   ! Correcting for 1/(microohm centimeter) units
  sigma = aux * sigma                                    ! Change sigma to 1/(microohm centimeter)
  L11   = aux * sigma                                    ! L11 is sigma   
  L12   = aux * L12                                      ! change L12 units in same way as sigma; now eV/uohm-cm
  L12   = L12 * 1.602176487E-19_DP                       ! fixing units  from eV/uohm-cm to Joule/uohm-cm
  thermopower1 = L12 / L11 / temperature/ELECTRON_SI     ! now Joule/(uohm-com) * (uohm-cm) / K / Coulomb = J/K/C = Volts/Kelvin
  L21   = aux * L21                                      ! change L21 units in same way as sigma; now eV/uohm-cm
  L21   = L21 * 1.602176487E-19_DP                       ! fixing units  from eV/uohm-cm to Joule/uohm-cm                   
  thermopower2 = L21 / L11 / temperature /ELECTRON_SI    ! now Joule/(uohm-com) * (uohm-cm) / K / Coulomb = J/K/C = Volts/Kelvin   
  thermopower  = (L21 + L12)/2.0_DP                      ! wont work like this....doing it by hand in the printout....
  L12_symmetric   = aux * L12_symmetric                  ! change L12_symmetric in same way as sigma; now eV/uohm-cm
  L12_symmetric   = L12_symmetric * 1.602176487E-19_DP   ! fixing units  from eV/uohm-cm to Joule/uohm-cm                  
  thermopower_symmetric = L12_symmetric/ L11 / temperature /ELECTRON_SI ! now Joule/(uohm-com) * (uohm-cm) / K / Coulomb = J/K/C = Volts/Kelvin
  L22aa   = aux * L22aa                                  ! change L22 in same way as sigma
  L22aa   = L22aa * 1.0E8_DP                             ! 1.0E8_DP is fixing units  from 1/uohm-cm to 1/ohm-m
  L22aa   = L22aa * 1.602176487E-19_DP                   ! change eV to Joule                      
  L22aa   = L22aa * 1.602176487E-19_DP                   ! change second eV to Joule   
  thermal_conductivity = (L22aa - (L12 * L21 / (sigma * 1.0E8_DP))) * 1/(ELECTRON_SI*ELECTRON_SI*temperature)

  ! performing simple integral for DC conductivity estimation
  conductivity_sum = 0.0_DP
  do iw = 1, nw
     conductivity_sum = conductivity_sum + (dw * fermi_distribution_derivative_neg(iw) * (real( (sigma(1,1,iw)+sigma(2,2,iw)+sigma(3,3,iw)) , DP )/3.0_DP))
  enddo

  ! Output

  open(11,file='optcond.dat', action='write')
  write(11,'(A)') '# Electric conductivity tensor in the xyz system of reference'
  write(11,*) '# Rough estimate of electrical conductivity in 1/(microohm*cm)', 2.0_DP * conductivity_sum
  write(11,'(A)'             ) '#      hbar*w (eV)                           sigma(w) ( 1/(microohm*cm) )'
  write(11,'(A)',advance='no') '# eV                               11                  12                  13'
  write(11,'(A)',advance='no') '                  21                  22                  23                  31'
  write(11,'(A)',advance='no') '                  32                  33               trace/3      '
  write(11,'(A)',advance='no') 'thermopowerL12(experimental)  thermopowerL21(experimental)  thermopowerL12_L21(experimental) '
  write(11,'(A)') 'thermopower_symmetric(experimental) thermal_conductivity(very_experimental)'

  do iw = 1, nw
     write(11,'(16(f40.20,1x))') wgrid(iw), ((real(sigma(ix1,ix2,iw),DP), ix2=1,3),ix1=1,3), &
                                        real( (sigma(1,1,iw)+sigma(2,2,iw)+sigma(3,3,iw)) , DP )/3.0_DP,&
                                        real( (thermopower1(1,1,iw)+thermopower1(2,2,iw)+thermopower1(3,3,iw)) , DP )/3.0_DP, &
                                        real( (thermopower2(1,1,iw)+thermopower2(2,2,iw)+thermopower2(3,3,iw)) , DP )/3.0_DP, &
                                        (real( (thermopower1(1,1,iw)+thermopower1(2,2,iw)+thermopower1(3,3,iw)) , DP )/3.0_DP +&
                                         real( (thermopower2(1,1,iw)+thermopower2(2,2,iw)+thermopower2(3,3,iw)) , DP )/3.0_DP)/2.0_DP, &
                                        real( (thermopower_symmetric(1,1,iw)+thermopower_symmetric(2,2,iw)+thermopower_symmetric(3,3,iw)) , DP )/3.0_DP, &
                                        real( (thermal_conductivity(1,1,iw)+thermal_conductivity(2,2,iw)+thermal_conductivity(3,3,iw)) , DP )/3.0_DP
  enddo 
  close(11)

!  open(11,file='optcond-xyz-imag.dat', action='write')
!  write(11,*)'# hbar*w (eV)   sigma(w) ( 1/(10^2 micro-ohm*cm) x 10^-6 )'
!  do iw = 1, nw
!     write(11,'(10f20.12)') wgrid(iw), ((aimag(sigma(ix1,ix2,iw)), ix2=1,3),ix1=1,3)
!  enddo
!  close(11)



  write(*,*) 'Done'

end program kubocalc
