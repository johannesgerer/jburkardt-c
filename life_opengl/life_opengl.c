/*
 * Conway's Game of Life in OpenGL
 * Simon Green (simes@hotbot.com) 3/00
 *
 * This uses roughly the same algorithm described in the OpenGL
 * red book (p.590 in 3rd edition) to compute Life in the frame buffer:
 *
 * - Reads frame buffer contents to texture using CopyTexImage.
 * - Renders textured quad 8 times offset one pixel in each direction,
 *   using alpha test and stencil buffer increment to count neighbours.
 * - Uses stencil test to determine deaths and births for next generation
 *
 * This is not very efficient (10 passes), but leaves results in texture map.
 * Requires stencil buffer, destination alpha.
 *
 */

# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>

/*
  For the MAC OSX, you need "GLUT/glut.h" rather than "GL/glut.h".
*/
# include <GLUT/glut.h>

int gen = 0;
int lifeh = 256;
int lifew = 256;
int winh = 256;
int winw = 256;

GLboolean animating = 0;
GLboolean wraparound = 1;

GLfloat alivecol[] = { 0.0, 1.0, 0.0, 1.0 };
GLfloat deadcol[] = { 0.0, 0.0, 0.0, 0.0 };

/*
   The R-Pentomino:
    OO
   OO
    O
*/

int rpentomino[] = {
    1, 0,
    0, 1,
    1, 1,
    1, 2,
    2, 2,
};

int main ( int argc, char **argv );
void checkErrors ( char *s );
void display ( void );
void drawPattern ( int *pattern, int n, int x, int y );
void drawQuad ( int x, int y );
void generation ( void );
void idle ( void );
void initGL ( void );
void key ( unsigned char key, int x, int y );
void reshape ( int w, int h );

/******************************************************************************/

int main ( int argc, char **argv )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LIFE_OPENGL.

  Author:

    Simon Green
*/
{
  printf ( "\n" );
  printf ( "LIFE_OPENGL:\n" );
  printf ( "  C version\n" );

    glutInit(&argc, argv);

    initGL();

    glutMainLoop();

    return 0;
}
/******************************************************************************/

void checkErrors ( char *s )

/******************************************************************************/
/*
  Purpose:

    CHECKERRORS ???
*/
{
  GLenum error;

  while ( ( error = glGetError ( ) ) != GL_NO_ERROR )
  {
    fprintf(stderr, "%s: error - %s\n", s, (char *) gluErrorString(error));
  }
  return;
}
/******************************************************************************/

void display ( void )

/******************************************************************************/
/*
  Purpose:

    DISPLAY ???.
*/
{
    generation();
    glFinish();

  return;
}
/******************************************************************************/

void drawPattern ( int *pattern, int n, int x, int y )

/******************************************************************************/
/*
  Purpose:

    DRAWPATTERN ???
*/
{
    int i;

    glPointSize(1.0);

    glPushMatrix();
    glLoadIdentity();
    glTranslatef((GLfloat) x, (GLfloat) y, 0.0);

    glColor4fv(alivecol);
    glBegin(GL_POINTS);

    for(i=0; i<n; i++)
    {
        glVertex2i(*pattern, *(pattern+1));
        pattern+=2;
    }
    glEnd();

    glPopMatrix();

  return;
}
/******************************************************************************/

void drawQuad ( int x, int y )

/******************************************************************************/
/*
  Purpose:

    DRAWQUAD ???
*/
{
    glBegin(GL_QUADS);

    if (wraparound) {
        // this is clever - use texture coordinate wrapping for edge wraparound
        float tx, ty;
        tx = x / (float) lifew;
        ty = y / (float) lifeh;
        glTexCoord2f(tx, ty);               glVertex2i(0, 0);
        glTexCoord2f(tx + 1.0, ty);         glVertex2i(lifew, 0);
        glTexCoord2f(tx + 1.0, ty + 1.0);   glVertex2i(lifew, lifeh);
        glTexCoord2f(tx, ty + 1.0);         glVertex2i(0, lifeh);

    } else {
        glTexCoord2f(0.0, 0.0); glVertex2i(x, y);
        glTexCoord2f(1.0, 0.0); glVertex2i(x + lifew, y);
        glTexCoord2f(1.0, 1.0); glVertex2i(x + lifew, y + lifeh);
        glTexCoord2f(0.0, 1.0); glVertex2i(x, y + lifeh);
    }

    glEnd();

  return;
}
/******************************************************************************/

void generation ( void )

/******************************************************************************/
/*
  Purpose:

    GENERATION computes the next generation.

  Discussion:

    There are 8 + 2 = 10 passes total.
*/
{
    glClear(GL_STENCIL_BUFFER_BIT);

/* 
  read frame buffer to texture memory 
*/
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, 0, 0, lifew, lifeh, 0);

    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
/* 
  Draw image 8 times offset with stencil to count neighbours 
*/
    glAlphaFunc(GL_GREATER, 0.0);   /* only write pixels with alpha > 0.0 */
    glEnable(GL_ALPHA_TEST);
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);    /* no need to write to color buffer */

    glStencilFunc(GL_ALWAYS, 0, 0xff);
    glStencilOp(GL_KEEP, GL_KEEP, GL_INCR); /* fail, zfail, zpass */
    glEnable(GL_STENCIL_TEST);

    drawQuad(-1, -1);
    drawQuad( 0, -1);
    drawQuad( 1, -1);
    drawQuad(-1,  0);
    drawQuad( 1,  0);
    drawQuad(-1,  1);
    drawQuad( 0,  1);
    drawQuad( 1,  1);

    glDisable(GL_STENCIL_TEST);
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_TEXTURE_2D);

    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);    /* enable writes to color buffer */

    // draw births (neighbours == 3)
    glStencilFunc(GL_EQUAL, 3, 0xff);       /* write and decr, only if stencil == 3 */
    glStencilOp(GL_KEEP, GL_DECR, GL_DECR); /* fail, zfail, zpass */
    glEnable(GL_STENCIL_TEST);

    glColor4fv(alivecol);
    drawQuad(0, 0);

    // clear everything but survivors (neighbours == 2 or 3) */
    glStencilFunc(GL_NOTEQUAL, 2, 0xff);    /* write only if stencil != 2 */
    glColor4fv(deadcol);
    drawQuad(0, 0);

    glDisable(GL_STENCIL_TEST);

    checkErrors("generation");
    gen++;

  return;
}

/******************************************************************************/

void idle ( void )

/******************************************************************************/
/*
  Purpose:

    IDLE ???.
*/
{
  glutPostRedisplay();

  return;
}
/******************************************************************************/

void initGL ( void )

/******************************************************************************/
/*
  Purpose:

    INITGL initializes OpenGL.
*/
{
  glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_STENCIL | GLUT_SINGLE);
  glutInitWindowSize(winw, winh);
  glutCreateWindow("OpenGL Game of Life");

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutIdleFunc(idle);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);

  glReadBuffer(GL_FRONT);
  glDrawBuffer(GL_FRONT);

  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClearStencil(0);
  glClear(GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
  reshape(lifew, lifeh);
  drawPattern(rpentomino, 5, lifew/2, lifeh/2);

  return;
}
/******************************************************************************/

void key ( unsigned char key, int x, int y )

/******************************************************************************/
/*
  Purpose:

    KEY ???.
*/
{
    switch (key) {
    case '\033':
        exit(0);
        break;

    case 'g':
        if (animating) {
            glutIdleFunc(NULL);
            animating = 0;
        } else {
            glutIdleFunc(idle);
            animating = 1;
        }
        break;

    case ' ':
        animating = 0;
        glutIdleFunc(NULL);
        glutPostRedisplay();
        break;
    }

  return;
}
/******************************************************************************/

void reshape ( int w, int h )

/******************************************************************************/
/*
  Purpose:

    RESHAPE ???.
*/
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, (GLdouble) w, 0.0, (GLdouble) h, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);

  glViewport(0.0f, 0.0f, w, h);

  winw = w;
  winh = h;

  return;
}


