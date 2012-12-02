
/* Two-Dimensional Sierpinski Gasket          */
/* Generated Using Randomly Selected Vertices */
/* And Bisection                              */

/* 
  #include <GL/glut.h>
*/
# include <GLUT/glut.h>

void myinit(void)
{
 
/* attributes */

      glClearColor(1.0, 1.0, 1.0, 1.0); /* white background */
      glColor3f(1.0, 0.0, 0.0); /* draw in red */

/* set up viewing */
/* 500 x 500 window with origin lower left */

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0.0, 500.0, 0.0, 500.0);
      glMatrixMode(GL_MODELVIEW);
}

void display( void )
{

/* define a point data type */

    typedef GLfloat point2[2];     

    point2 vertices[3]={{0.0,0.0},{250.0,500.0},{500.0,0.0}}; /* A triangle */

    int i, j, k;
    int rand();       /* standard random number generator */
    point2 p ={75.0,50.0};  /* An arbitrary initial point inside traingle */

    glClear(GL_COLOR_BUFFER_BIT);  /*clear the window */


/* compute and plots 50000 new points */

    for( k=0; k<50000; k++)
    {
         j=rand()%3; /* pick a vertex at random */


     /* Compute point halfway between selected vertex and old point */

         p[0] = (p[0]+vertices[j][0])/2.0; 
         p[1] = (p[1]+vertices[j][1])/2.0;
   
     /* plot new point */

          glBegin(GL_POINTS);
               glVertex2fv(p); 
          glEnd();
   
     }
     glFlush(); /* clear buffers */
 }

void main(int argc, char** argv)
{

/* Standard GLUT initialization */

    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB); /* default, not needed */
    glutInitWindowSize(500,500); /* 500 x 500 pixel window */
    glutInitWindowPosition(0,0); /* place window top left on display */
    glutCreateWindow("Sierpinski Gasket"); /* window title */
    glutDisplayFunc(display); /* display callback invoked when window opened */

    myinit(); /* set attributes */

    glutMainLoop(); /* enter event loop */
}

