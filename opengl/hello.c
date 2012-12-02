# include <stdlib.h>

# include <GLUT/glut.h>
//#include <GL/glut.h>

int main ( int argc, char** argv );
void display ( void );
void init ( void );

/******************************************************************************/

int main ( int argc, char** argv )

/******************************************************************************/
{

/* 
  Pass command line arguments to GLUT.
*/
  glutInit ( &argc, argv );

/* 
  Declare display mode (single buffer and RGBA).
*/
  glutInitDisplayMode ( GLUT_SINGLE | GLUT_RGB );

/* 
  Declare initial window size and position.
*/
  glutInitWindowSize ( 250, 250 ); 
  glutInitWindowPosition ( 100, 100 );

/* 
  Open window with "hello" in its title bar.  
*/
  glutCreateWindow ( "hello" );

/* 
  Call initialization routines.
*/
  init ();

/* 
  Register callback function to display graphics.
*/
  glutDisplayFunc ( display ); 

/* 
  Enter main loop and process events.
*/
  glutMainLoop();

  return 0;
}
/******************************************************************************/

void display ( void )

/******************************************************************************/
{

/* 
  Clear all pixels.
*/
  glClear ( GL_COLOR_BUFFER_BIT );
/* 
  Set the color.  
*/
  glColor3f ( 0.1, 0.3, 1.0 );
/* 
  Draw the polygon with corners at
 (0.25, 0.25, 0.0) and (0.75, 0.75, 0.0)  
*/
  glBegin ( GL_POLYGON );
    glVertex3f ( 0.25, 0.25, 0.0 );
    glVertex3f ( 0.75, 0.25, 0.0 );
    glVertex3f ( 0.75, 0.75, 0.0 );
    glVertex3f ( 0.25, 0.75, 0.0 );
  glEnd();

/* 
  Don't wait!  
  start processing buffered OpenGL routines.
*/
  glFlush ();

  return;
}
/******************************************************************************/

void init ( void ) 

/******************************************************************************/
{
/* 
  Select the clearing color.
*/
  glClearColor ( 0.0, 0.0, 0.0, 0.0 );

/* 
  Initialize viewing values.
*/
  glMatrixMode ( GL_PROJECTION );

  glLoadIdentity();

  glOrtho ( 0.0, 1.0, 0.0, 1.0, -1.0, 1.0 );

  return;
}


