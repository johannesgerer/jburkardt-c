c  "Floating comments" are assigned to the next module.
c
      program alpha
c
c*********************************************************************72
c
cc ALPHA is a main program.
c
      call beta ( y )
      x = gamma ( y )
      stop
      end
c
c*********************************************************************72
c
c  "Nameless" code like this is assumed to be a (main) program.
c
      write ( *, * ) 'This is a main program.'

      stop
      end
      subroutine beta ( y )
c
c*********************************************************************72
c
cc BETA is a subroutine
c
      y = 3.14159265

      return
      end
      function gamma ( x )
c
c*********************************************************************72
c
cc GAMMA is a function with no type statement.
c
      real gamma
      real x

      gamma = sqrt ( x )

      return
      end
      complex function delta ( x )
c
c*********************************************************************72
c
cc DELTA is a complex function with type statement.
c
      complex x

      delta = x * x

      return
      end
      integer function epsilon ( x )
c
c*********************************************************************72
c
cc EPSILON is an integer function with type statement.
c
      integer x

      epsilon = x + 1

      return
      end
c
c  Inter-module comments go to the next module.
c
      real function zeta ( x )
c
c*********************************************************************72
c
cc ZETA is a real function with type statement.
c
      real x

      zeta = 1.0 / x

      return
      end
      module eta
c
c*********************************************************************72
c
cc ETA is a module.
c
      real x

      x = 17.0

      end
      block data theta
c
c*********************************************************************72
c
cc THETA is a block data routine.
c
      real x
      common / samantha / x
      save / samantha /
      data x / 1.0 /

      end
      procedure sigma ( w )
c
c*********************************************************************72
c
cc SIGMA is a procedure, which is not a kosher FORTRAN module type.
c
      real w

      w = w * 2

      end
c
c  Trailing comments go where, exactly?
c
