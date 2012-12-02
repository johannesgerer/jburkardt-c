#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vec.h"
#include "boid.h"

#define rrand(a) (rand()%a)

#define DEFAULT_CENTER_BIAS 7
#define DEFAULT_AVG_VEL     3
#define DEFAULT_CHILLING    1

Boid new_boid(int W, int H)
{
   Boid boid;
   double px,py,pz;
   double vx,vy,vz;
   
   boid=(Boid)malloc(sizeof(_Boid));
   
   px=(double)(rrand((W<<4)) - (W<<3));
   py=(double)(rrand((H<<4)) - (H<<3));
   pz=(double)(rrand(2000) + 2000);
   
   boid->pos=new_vec(px,py,pz);
   
   vx=(double)(rrand(51) - 25);
   vy=(double)(rrand(51) - 25);
   vz=(double)(rrand(51) - 25);
   
   boid->vel=new_vec(vx,vy,vz);
   
   boid->wing_level=(int)(rrand(200)-100);
   
   boid_perspective(boid, W, H);
   
   return boid;
}

void boid_perspective(Boid boid, int W, int H)
{
   double zfactor, zf;
   Vec tail, tail_end;
   double tail_lx, tail_lz, tail_rx, tail_rz;
   double tailx, tailz;
         
   tail=vec_copy(boid->vel);
   tail_end=vec_copy(boid->vel);
   
   if(boid->pos->z <= 0) {
      boid->onscreen = 0;
   }
   else {
      zf = W/(double)2.5;
      zfactor=((double)boid->pos->z)/zf;
      
      boid->X = (W>>1) + (int)(boid->pos->x/zfactor);
      boid->Y = (H>>1) + (int)(boid->pos->y/zfactor);
            
      boid->shadow[0].x = boid->X;
      boid->shadow[0].y = (H>>1) + (int)(1000/zfactor);
      
      vec_setmag(tail_end, 40);
      vec_diff(boid->pos, tail_end, tail_end);
      
      zfactor=((double)tail_end->z)/zf;
      boid->tail_X = (W>>1) + (int)(tail_end->x/zfactor);
      boid->tail_Y = (H>>1) + (int)(tail_end->y/zfactor);
      boid->shadow[2].x = boid->tail_X;
      boid->shadow[2].y = (H>>1) + (int)(1000/zfactor);
      
      vec_setmag(tail, 50);
      vec_diff(boid->pos, tail, tail);
      
      tailx = -tail->z/60;
      tailz = tail->x/60;
      
      tail_lx = tail->x - tailx;
      tail_lz = tail->z - tailz;
      
      tail_rx = tail->x + tailx;
      tail_rz = tail->z + tailz;
      
      tail->y -= boid->wing_level;
      
      zfactor = ((double)tail_lz)/zf;
      boid->tail_lX = (W>>1) + (int)(tail_lx/zfactor);
      boid->tail_lY = (H>>1) + (int)(tail->y/zfactor);
      boid->shadow[1].x = boid->tail_lX;
      boid->shadow[1].y = (H>>1) + (int)(1000/zfactor);
      
      zfactor = ((double)tail_rz)/zf;
      boid->tail_rX = (W>>1) + (int)(tail_rx/zfactor);
      boid->tail_rY = (H>>1) + (int)(tail->y/zfactor);
      boid->shadow[3].x = boid->tail_rX;
      boid->shadow[3].y = (H>>1) + (int)(1000/zfactor); 
      
      
      boid->onscreen = boid_isonscreen(boid,W,H);
   }
   
   free(tail);
   free(tail_end);
}

int boid_isonscreen(Boid boid, int W, int H)
{
   return (boid->X>=0 && boid->X<W && boid->Y>=0 && boid->Y<H);
}

/*
 * Move this boid, wrt allboids
 */
void boid_move(Boid boid, Boid allboids[], int numboids,
	       Vec real_center, Vec real_avgvel,
	       int W, int H)
{
   Vec center, center_bias;
   Vec avgvelocity, avgvel_bias;
   Vec chilling;
   
   if(boid->perching) {
      if(boid->perch_timer > 0) {
	 boid->perch_timer--;
	 return;
      }
      else {
	 boid->perching = 0;
      }
   }

   center_bias = zero_vec();
   center = boid_perceive_center(boid, real_center, numboids);
   vec_diff(center, boid->pos, center_bias);
   vec_rshift(center_bias, DEFAULT_CENTER_BIAS);
   
   avgvel_bias = zero_vec();
   avgvelocity = boid_av_vel(boid, real_avgvel, numboids);
   vec_diff(avgvelocity, boid->vel, avgvel_bias);
   vec_rshift(avgvel_bias, DEFAULT_AVG_VEL);
   
   chilling = boid_chill_out(boid, allboids, numboids);
   vec_rshift(chilling, DEFAULT_CHILLING);
   
   vec_add(boid->vel, center_bias);
   vec_add(boid->vel, avgvel_bias);
   vec_add(boid->vel, chilling);
   
   vec_limit(boid->vel, 100);
   
   vec_add(boid->pos, boid->vel);

   free(center); free(center_bias);
   free(avgvelocity); free(avgvel_bias);
   free(chilling);
   
   if(boid->upstroke) {
      if(boid->wing_level >= 100) {	  
	 boid->upstroke = 0;
      }
      else {
	 boid->wing_level += 40;
      }
   }
   else {
      if(boid->wing_level <= -100) {
	 boid->upstroke = 1;
      }
      else {
	 boid->wing_level -= 20;
      }
      
   }
   
   
   /* bound world */
   if(boid->pos->x < -1500) {
      boid->vel->x += 10;
   }
   if(boid->pos->x > 1500) {
      boid->vel->x -= 10;
   }
   if(boid->pos->y < -1200) {
      boid->vel->y += 10;
   }
   if(boid->pos->y > 800) {
      boid->vel->y -= 10;
   }
   if(boid->pos->y > 1000) {
      boid->pos->y = 1000; /* Hit ground!! */
      boid->perching = 1;
      boid->perch_timer = (int)(rrand(20) + 30);
      boid->wing_level = 60;
      boid->vel->y = 0;
   }
   if(boid->pos->z < 1500) {
      boid->vel->z += 10;
   }
   if(boid->pos->z > 3000) {
      boid->vel->z -= 10;
   }
      
   boid_perspective(boid, W, H);
}

Vec boid_perceive_center(Boid boid, Vec real_cent, int numboids)
{
   Vec perc_cent;
   
   perc_cent = zero_vec();
   vec_diff(real_cent, boid->pos, perc_cent);
   vec_sdiv(perc_cent, (double)(numboids-1));
   
   return perc_cent;
}

Vec boid_av_vel(Boid boid, Vec real_avgvel, int numboids)
{
   Vec perc_avgvel;
   
   perc_avgvel = zero_vec();
   vec_diff(real_avgvel, boid->vel, perc_avgvel);
   vec_sdiv(perc_avgvel, (double)(numboids-1));
   
   return perc_avgvel;
}

Vec boid_chill_out(Boid boid, Boid boids[], int numboids)
{
   Vec chill, bigchill;
   int i;
   
   chill = zero_vec();
   bigchill = zero_vec();
   
   for(i=0; i<numboids; i++) {
      if(boids[i] != boid) {
	 if(vec_rdist(boid->pos, boids[i]->pos) <= 100) {
	    vec_diff(boid->pos, boids[i]->pos, chill);
	    vec_add(bigchill, chill);
	 }
      }
   }
   
   free(chill);
   
   return bigchill;
}

