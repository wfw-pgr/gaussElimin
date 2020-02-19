! ====================================================== !
! === Gauss Elimination Example                      === !
! ====================================================== !
program main
  use gaussElimMod, only : gaussElimin
  implicit none
  integer         , parameter   :: nSize    = 3
  character(100)  , parameter   :: FileName = 'dat/eq.inp'
  double precision, allocatable :: Amat(:,:), bvec(:), xvec(:), res(:)
  integer                       :: i, j
  double precision              :: avg, std

  ! ------------------------------------- !
  ! --- [1] Data Load                 --- !
  ! ------------------------------------- !
  allocate( Amat(nSize,nSize), bvec(nSize), xvec(nSize) )
  open(50,file=trim(FileName),status='old',form='formatted')
  do i=1, nSize
     read(50,*) Amat(i,:), bvec(i)
  enddo
  close(50)

  do i=1, nSize
     write(6,*) Amat(i,:), '|', bvec(i)
  enddo

  ! ------------------------------------- !
  ! --- [2] Gauss Elimination         --- !
  ! ------------------------------------- !
  call gaussElimin( Amat, xvec, bvec, nSize )

  write(6,*) xvec(:)
  
  ! ------------------------------------- !
  ! --- [3] Answer Check              --- !
  ! ------------------------------------- !
  allocate( res(nSize) )
  do i=1, nSize
     res(i) = 0.d0
     do j=1, nSize
        res(i) = res(i) + Amat(i,j) * xvec(j)
     enddo
  enddo
  res(:) = bvec(:) - res(:)
  avg = 0.d0
  std = 0.d0
  do i=1, nSize
     avg = avg + res(i)
     std = std + res(i)**2
     write(6,*) res(i)
  enddo
  avg = avg / dble( nSize )
  std = std / dble( nSize )
  std = sqrt( std - avg**2 )
  write(6,*) 'avg :: ', avg, 'std :: ', std
  
  return
end program main
