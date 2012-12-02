!  "Floating comments" are assigned to the next module.
!
program alpha
!
!*****************************************************************************80
!
!! ALPHA is a main program.
!
  call beta ( y )
  x = gamma ( y )
  stop
end
!
!*****************************************************************************80
!
!  "Nameless" code like this is assumed to be a (main) program.
!
  write ( *, * ) 'This is a main program.'

  stop
end
subroutine beta ( y )
!
!*****************************************************************************80
!
!! BETA is a subroutine
!
  y = 3.14159265

  return
end
function gamma ( x )
!
!*****************************************************************************80
!
!! GAMMA is a function with no type statement.
!
  real gamma
  real x

  gamma = sqrt ( x )

  return
end
complex function delta ( x )
!
!*****************************************************************************80
!
!! DELTA is a complex function with type statement.
!
  complex x

  delta = x * x

  return
end
integer function epsilon ( x )
!
!*****************************************************************************80
!
!! EPSILON is an integer function with type statement.
!
  integer x

  epsilon = x + 1

  return
end
!
!  Inter-module comments go to the next module.
!
real function zeta ( x )
!
!*****************************************************************************80
!
!! ZETA is a real function with type statement.
!
  real x

  zeta = 1.0 / x

  return
end
module eta
!
!*****************************************************************************80
!
!! ETA is a module.
!
  real x

  x = 17.0

end
block data theta
!
!*****************************************************************************80
!
!! THETA is a block data routine.
!
  real x
  common / samantha / x
  save / samantha /
  data x / 1.0 /

end
procedure sigma ( w )
!
!*****************************************************************************80
!
!! SIGMA is a procedure, which is not a kosher FORTRAN module type.
!
  real w

  w = w * 2

end
!
!  Trailing comments go where, exactly?
!
