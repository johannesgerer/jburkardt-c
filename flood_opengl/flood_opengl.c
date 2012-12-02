/* Asn2 flooding Truchet tiles IS2780 computer graphics 9/25/2001 A.Wetzel */
/* Goal1: understand that images are just numerical patterns in memory */
/* Goal2: show transfer of data from CPU space to the graphics system */
/* Goal3: introduce other basic OpenGL operations and mouse handling */

# include <GLUT/glut.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# define PI 3.14159265358979323846

void glut_setup(void) ;
void gl_setup(void) ;
void display(void) ;
void display2(void) ;
void idle(void) ;
void keyboard(unsigned char key, int mouse_x, int mouse_y) ;
void mouse( int button, int state, int x, int y) ;
void data_setup(void) ;
void draw_tile_a(void) ;
void draw_tile_b(void) ;
void arc(float r, float start_deg, float end_deg) ;
void putpixel(int x, int y, int r, int g, int b) ;
void getpixel(int x, int y, int *r, int *g, int *b) ;
void draw_pattern(void) ;
void do_floodfill(int x, int y) ;
void floodfill(int x, int y, int r, int g, int b);

/* 
  opengl/glut info 
*/
int window_x_size = 512;
int window_y_size = 512;
int main_window_id;
int second_window_id;
int first_draw_flag = 0;

# define MAX_ROW 32
# define MAX_COL 32

int tile_pattern[MAX_ROW][MAX_COL];

# define PIXEL_ROW 512
# define PIXEL_COL 512

char pixels[PIXEL_ROW][PIXEL_COL][3];

/*******************************************************************************/

int main ( int argc, char **argv ) 

/*******************************************************************************/
/*
  Purpose:

    MAIN is the main program for FLOOD_OPENGL.

  Author:

    Art Wetzel
*/
{
  glutInit ( &argc, argv );

  glut_setup ( );

  gl_setup ( );

  data_setup ( );

  glutMainLoop ( );

  return 0;
}
/*******************************************************************************/

void glut_setup ( void ) 

/*******************************************************************************/
/*
  Purpose:

    GLUT_SETUP

  Author:

    Art Wetzel
*/
{
  glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA);
  glutInitWindowSize(window_x_size, window_y_size);

  main_window_id = glutCreateWindow("Main Window");
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutIdleFunc(idle);
  gl_setup();

  second_window_id = glutCreateWindow("Second Window");
  glutDisplayFunc(display2);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  glutMouseFunc(mouse);
  gl_setup();

  glutSetWindow(main_window_id);

  return;
}
/*******************************************************************************/

void gl_setup ( void ) 

/*******************************************************************************/
/*
  Purpose:

    GL_SETUP

  Author:

    Art Wetzel
*/
{
  glClearColor ( 1, 1, 1, 0 );

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, window_x_size, 0, window_y_size);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt ( 0, 0, 1, 0, 0, 0, 0, 1, 0 );

  return;
}
/*******************************************************************************/

void display ( void ) 

/*******************************************************************************/
/*
  Purpose:

    DISPLAY generates the graphics output.

  Author:

    Art Wetzel
*/
{
  glClear(GL_COLOR_BUFFER_BIT);

  glPushMatrix();
  glScalef(16, 16, 1);
  draw_pattern();
  glPopMatrix();

  glutSwapBuffers();

  if(first_draw_flag == 0) 
  {
    first_draw_flag = 1;
    glReadPixels(0, 0, 512, 512, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  }
  return;
}
/*******************************************************************************/

void display2 ( void ) 

/*******************************************************************************/
/*
  Purpose:

    DISPLAY2

  Author:

    Art Wetzel
*/
{
  glClear(GL_COLOR_BUFFER_BIT);
  glRasterPos2i(0, 0);
  glDrawPixels(512, 512, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  glutSwapBuffers();

  return;
}
/*******************************************************************************/

void idle ( void ) 

/*******************************************************************************/
/*
  Purpose:

    IDLE

  Author:

    Art Wetzel
*/
{
  glutSetWindow(main_window_id);
  glutPostRedisplay();

  glutSetWindow(second_window_id);
  glutPostRedisplay();

  return;
}
/*******************************************************************************/

void keyboard ( unsigned char key, int mouse_x, int mouse_y ) 

/*******************************************************************************/
/*
  Purpose:

    KEYBOARD reacts to keyboard events.

  Author:

    Art Wetzel
*/
{
  switch(key)
  {
    case 27:
      exit(0);
      break;

    default:
      printf("got key %c, mouse at %d %d\n",key, mouse_x, mouse_y);
      break;
  }

  return;
}
/*******************************************************************************/

void mouse ( int button, int state, int x, int y ) 

/*******************************************************************************/
/*
  Purpose:

    MOUSE reacts to mouse events.

  Author:

    Art Wetzel
*/
{
  printf("got mouse %d @ %d %d\n",button,x,y);

/* 
  map mouse y to array y 
*/
  y = 512 - y;

  switch(button) 
  {
    case GLUT_LEFT_BUTTON:
      if ( state == GLUT_DOWN ) 
      {
        do_floodfill ( x, y );
      }
      break;
  }
  return;
}
/*******************************************************************************/

void data_setup ( void ) 

/*******************************************************************************/
/*
  Purpose:

    DATA_SETUP

  Author:

    Art Wetzel
*/
{
  int i, j;

  for(i = 0; i < MAX_ROW; i++)
  {
    for(j = 0; j < MAX_COL; j++)
    {
      tile_pattern[i][j] = rand() & 1;
    }
  }
  return;
}
/*******************************************************************************/

void draw_tile_a ( void ) 

/*******************************************************************************/
/*
  Purpose:

    DRAW_TILE_A

  Author:

    Art Wetzel
*/
{
  /*
  * hint: start by drawing an empty box
  * and positioning it correctly
  */
  glColor3f(0, 0, 0);

  glPushMatrix();
  glTranslatef(0, 1, 0);
  arc(0.5, 270.0, 360.0);
  glPopMatrix();

  glPushMatrix();
  glTranslatef(1, 0, 0);
  arc(0.5, 90.0, 180.0);
  glPopMatrix();

  return;
}
/*******************************************************************************/

void draw_tile_b ( void ) 

/*******************************************************************************/
/*
  Purpose:

    DRAW_TILE_B

  Author:

    Art Wetzel
*/
{
  glColor3f(0, 0, 0);

  glPushMatrix();
  glTranslatef(0, 0, 0);
  arc(0.5, 0.0, 90.0);
  glPopMatrix();

  glPushMatrix();
  glTranslatef(1, 1, 0);
  arc(0.5, 180.0, 270.0);
  glPopMatrix();

  return;
}
/*******************************************************************************/

void draw_pattern ( void )

/*******************************************************************************/
/*
  Purpose:

    DRAW_PATTERN

  Author:

    Art Wetzel
*/
{
  int i, j;

  for(i = 0; i < MAX_ROW; i++)
  {
    for(j = 0; j < MAX_COL; j++)
    {
      switch(tile_pattern[i][j])
      {
      case 0:
        draw_tile_a();
        break;

      case 1:
        draw_tile_b();
        break;

      default:
        printf("got   tile_pattern[i][j] %d\n",
          tile_pattern[i][j]);
        break;
      }
      glTranslatef(1, 0, 0);
    }
    glTranslatef(-MAX_COL, 0, 0);
    glTranslatef(0, 1, 0);
  }
  return;
}
/*******************************************************************************/

void arc ( float r, float start_deg, float end_deg )

/*******************************************************************************/
/*
  Purpose:

    ARC

  Author:

    Art Wetzel
*/
{
  float ang;

  glBegin ( GL_POINTS );

  for ( ang = start_deg * PI / 180.0; ang < end_deg * PI / 180.0; ang = ang + PI / 180.0 )
  {
    glVertex2f ( r * cos ( ang ), r * sin ( ang ) );
  }
  glEnd ( );

  return;
}
/*******************************************************************************/

void putpixel ( int x, int y, int r, int g, int b ) 

/*******************************************************************************/
/*
  Purpose:

    PUTPIXEL

  Author:

    Art Wetzel

  Parameters:

    Input, int X, Y, the row and column of the pixel.

*/
{
  char *p;
/*
  If the pixel is outside the range, return.
*/
  if ( ( x < 0 ) || 
       ( y < 0 ) || 
       ( PIXEL_COL <= x ) || 
       ( PIXEL_ROW <= y ) ) 
    return;
  }
/* 
  illustrate 2D addressing of 1D memory 
*/
  p = (char *)pixels + 3*x + y*(3*PIXEL_COL);

  *p = r; 
   p++;

  *p = g; 
   p++;

  *p = b;

  return;
}
/*******************************************************************************/

void getpixel ( int x, int y, int *r, int *g, int *b ) 

/*******************************************************************************/
/*
  Purpose:

    GETPIXEL returns the color of a pixel.

  Author:

    Art Wetzel

  Parameters:

    Input, int X, Y, the row and column of the pixel.

    Output, int *R, *G, *B, the current color of the pixel.
*/
{
  char *p;
/*
  If the pixel is outside the range, return black.
*/
  if ( ( x < 0 ) || 
       ( y < 0 ) || 
       ( PIXEL_COL <= x ) || 
       ( PIXEL_ROW <= y ) ) 
  {
    *r = 0;
    *g = 0;
    *b = 0;
  }
  else
  {
    p = ( char * ) pixels + 3 * x + y * ( 3 * PIXEL_ROW );
    *r = (unsigned char) *p++;
    *g = (unsigned char) *p++;
    *b = (unsigned char) *p;
  }

  return;
}
/*******************************************************************************/

void do_floodfill ( int x, int y ) 

/*******************************************************************************/
/*
  Purpose:

    DO_FLOODFILL chooses a color, and fills the selected region.

  Author:

    Art Wetzel

  Parameters:

    Input, int X, Y, the row and column of the pixel.
*/
{
  int r, g, b;
/*
  Choose a random 8 bit color.
*/
  r = rand ( ) & 0xFF;
  g = rand ( ) & 0xFF;
  b = rand ( ) & 0xFF;
/*
  Apply that color to this pixel, and all its neighbors.
*/
  floodfill ( x, y, r, g, b );

  return;
}
/*******************************************************************************/

void floodfill ( int x, int y, int r, int g, int b ) 

/*******************************************************************************/
/*
  Purpose:

    FLOODFILL recursively colors a pixel and all its neighbors.

  Author:

    Art Wetzel

  Parameters:

    Input, int X, Y, the row and column of the pixel.

    Input, int R, G, B, the color to be applied to the pixel and its neighbors.
*/
{
  int r1, g1, b1;
/*
  Retrieve the current color of the pixel at (X,Y).
*/
  getpixel ( x, y, &r1, &g1, &b1 );
/*
  Is the pixel white?
  If so, color it, and check its four immediate neighbors.
*/
  if ( r1 == 255 && g1 == 255 && b1 == 255 ) 
  {
    putpixel ( x, y, r, g, b );
    floodfill ( x+1, y, r, g, b );
    floodfill ( x-1, y, r, g, b );
    floodfill ( x, y+1, r, g, b );
    floodfill ( x, y-1, r, g, b );
  }
  return;
}
