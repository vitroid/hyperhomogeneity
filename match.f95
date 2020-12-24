subroutine match(qS, nS, qL, nL, thres2, err)
  implicit none
  ! match small and large point sets
  real(kind=8), intent(in)  :: qS(nS,3), qL(nL,3), thres2
  integer,      intent(in)  :: nS, nL
  real(kind=8), intent(out) :: err
  real(kind=8) :: d(3), dd
  integer      :: i, j
  logical      :: found

  err = 0
  do i=1, size(qS, 1)
     found = .FALSE.
     do j=1, size(qL, 1)
        d(:) = qS(i,:) - qL(j,:)
        dd   = d(1)**2 + d(2)**2 + d(3)**2
        if ( dd < thres2 ) then
           err = err + dd
           found = .TRUE.
           exit
        endif
     enddo
     if ( .not. found ) then
        err = -1
        return
     endif
  enddo
end subroutine match

