# include <stdlib.h>
# include <stdio.h>
# include <string.h>

int main ( );
void fred ( char **name );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for CHARACTER_ARG.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 May 2013

  Author:

    John Burkardt
*/
{
  char *name;

  printf ( "\n" );
  printf ( "CHARACTER_ARG:\n" );
  printf ( "  C version\n" );
  printf ( "  Demonstrate how a C function can return character data\n" );
  printf ( "  through the argument list.\n" );
  printf ( "\n" );
  printf ( "  Our main program declares a character pointer:\n" );
  printf ( "    char *name;\n" );
  printf ( "  then calls function fred() with the ADDRESS of the pointer:\n" );
  printf ( "    fred ( &name );\n" );
  printf ( "  Function fred receives its argument as\n" );
  printf ( "    void fred ( char **name )\n" );
  printf ( "  It allocates memory for *name:\n" );
  printf ( "    *name = ( char * ) malloc ( 12 * sizeof ( char ) );\n" );
  printf ( "  It sets a value to *name:\n" );
  printf ( "    strcpy ( *name, \"ob_data.txt\" );\n" );
  printf ( "  The main program now has a string stored in name:\n" );
  printf ( "    printf ( \"String returned from fred: %s\\n\", name )\n" );

  fred ( &name );

  printf ( "\n" );
  printf ( "  The value of name is now = \"%s\".\n", name );

  printf ( "\n" );
  printf ( "CHARACTER_ARG:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void fred ( char **name )

/******************************************************************************/
/*
  Purpose:

    FRED returns character data through one of its arguments.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 May 2013

  Author:

    John Burkardt

  Parameter:

    Input, char **NAME, the address of the address of a character.
    The value *NAME is the address of a character, and can be set
    to the address of a string of interest to the user.
*/
{
/*
  Create space, and store its address in *NAME.
*/
  *name = ( char * ) malloc ( 12 * sizeof ( char ) );
/*
  Park your data in the space.
*/
  strcpy ( *name, "ob_data.txt" );

  return;
}

