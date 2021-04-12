#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <math.h>
#include <iostream>
#include <chrono>
#include <openacc.h>

#include "Elastic.h"
#include "Plastic.h"

// set window dimensions
const static int WINDOW_WIDTH = 800;
const static int WINDOW_HEIGHT = 600;
const static double VIEW_WIDTH = 1.5*800.f;
const static double VIEW_HEIGHT = 1.5*600.f;
const static float H = 16.f; // kernel radius

// set time-stepping parameters
const float t_final = 10.0;
const int n_time_steps = 1000;
const float dt = t_final / float(n_time_steps);
const float idt = 1.0 / dt;

// define gravitation
const float gx = +0.0;
const float gy = -0.01;

// define material parameters
const float rho = 1.0;
const float E = 1000.0;
const float nu = 0.3;
const float yield = 20.0;

const float mu2 = E/(1.0+nu);
const float lam = mu2*nu/(1.0-2.0*nu);
const float kappa = lam+0.5f*mu2;
const float pmod = lam+mu2;
const float c = sqrt(pmod/rho);

// define material model
//Elastic* model = new Elastic(rho,E,nu);
//#pragma acc update device(model)
//Material* model = new Plastic(rho,E,nu,yield);
Plastic* model = new Plastic(rho,E,nu,yield);
//Elastic model = Elastic(rho,E,nu);

// define grid dimensions
const float dx = 10.0;
const float dy = 10.0;
const int n_elems_x = 50;
const int n_elems_y = 100;
const int n_elems = n_elems_x*n_elems_y;
const int n_nodes_x = n_elems_x + 1;
const int n_nodes_y = n_elems_y + 1;
const int n_nodes = n_nodes_x*n_nodes_y;

// define all nodes in the grid
static float x_0[n_nodes];
static float y_0[n_nodes];
static float x[n_nodes];
static float y[n_nodes];
static float m[n_nodes];
static float vx[n_nodes];
static float vy[n_nodes];
static float fx[n_nodes];
static float fy[n_nodes];

// define all elements in the grid
const int n_state_vars = model->numStateVariables();
static int connectivity[4*n_elems];
const int tot_state_vars = 4*n_elems*n_state_vars;
float* state_vars = new float[tot_state_vars];
static float plot[n_elems];
const float max_plot = 0.5;

// define element constants
const float XIa[4]  = { -1.0, +1.0, +1.0, -1.0 };
const float ETAa[4] = { -1.0, -1.0, +1.0, +1.0 };
const float gauss = 1.0 / sqrt(3.0);

void InitFEM(void) {

  std::cout << "setting up problem..." << std::endl;

  // copy material model pointer to device
  //#pragma acc enter data copyin(model[0:1])
  //acc_attach((void**)&model);

  std::cout << "defining initial nodal coordinates..." << std::endl;

  // define all nodes in the grid
  for (int iy = 0; iy < n_nodes_y; iy++) {
    for (int ix = 0; ix < n_nodes_x; ix++) {
      int id = n_nodes_x*iy+ix;
      x_0[id] = VIEW_WIDTH*0.3 + ix*dx;
      y_0[id] = VIEW_HEIGHT*0.15 + iy*dy;
      x[id] = x_0[id];
      y[id] = y_0[id];
      m[id] = 0.0;
      vx[id] = 0.0;
      vy[id] = 0.0;
      fx[id] = 0.0;
      fy[id] = 0.0;
    } // for ix = ...
  } // for iy = ...

  std::cout << "defining element connectivities..." << std::endl;

  // define all elements in the grid
  for (int iy = 0; iy < n_elems_y; iy++) {
    for (int ix = 0; ix < n_elems_x; ix++) {
      int id = 4*(n_elems_x*iy+ix);
      connectivity[id  ] = n_nodes_x*iy+ix;
      connectivity[id+1] = n_nodes_x*iy+ix+1;
      connectivity[id+2] = n_nodes_x*(iy+1)+ix+1;
      connectivity[id+3] = n_nodes_x*(iy+1)+ix;
      plot[n_elems_x*iy+ix] = 0.0;
      for (int iq = 0; iq < 4; iq++) {
	model->initialize(&state_vars[n_state_vars*(id+iq)]);
      } // for iq = ...
    } // for ix = ...
  } // for iy = ...

  std::cout << "computing lumped nodal masses..." << std::endl;

  // sum nodal masses
  const float dm = 0.25*rho*dx*dy;
  for (int ie = 0; ie < n_elems; ie++) {
    m[connectivity[4*ie  ]] += dm;
    m[connectivity[4*ie+1]] += dm;
    m[connectivity[4*ie+2]] += dm;
    m[connectivity[4*ie+3]] += dm;
  } // for ie = ...

  // invert nodal masses
  for (int in = 0; in < n_nodes; in++) {
    m[in] = 1.0 / m[in];
  } // for in = ...

} // InitFEM()

void TimeIntegrate(void) {

  float t = 0.0;

  std::cout << "starting loop over time steps..." << std::endl;

  auto t1 = std::chrono::high_resolution_clock::now();

#pragma acc data copy(x,y,vx,vy,plot,state_vars[0:tot_state_vars]) copyin(x_0,y_0,m,fx,fy,connectivity,m,gauss,XIa,ETAa,gx,gy,dt,idt,model)
  {
  
  // loop over all time steps
  while (t < t_final) {
    // update the current analysis time
    t += dt;

    // loop over all elements and compute internal forces
#pragma acc parallel loop async(1) present(x,y,vx,vy,plot,x_0,y_0,fx,fy,connectivity,state_vars) firstprivate(gauss,XIa,ETAa,idt)
    for (int ie = 0; ie < n_elems; ie++) {
	    
      // load local element coordinates from memory
      float xe[4] = { x[connectivity[4*ie  ]],
		      x[connectivity[4*ie+1]],
		      x[connectivity[4*ie+2]],
		      x[connectivity[4*ie+3]] };
      float ye[4] = { y[connectivity[4*ie  ]],
		      y[connectivity[4*ie+1]],
		      y[connectivity[4*ie+2]],
		      y[connectivity[4*ie+3]] };
      float x0e[4] = { x_0[connectivity[4*ie  ]],
		       x_0[connectivity[4*ie+1]],
		       x_0[connectivity[4*ie+2]],
		       x_0[connectivity[4*ie+3]] };
      float y0e[4] = { y_0[connectivity[4*ie  ]],
		       y_0[connectivity[4*ie+1]],
		       y_0[connectivity[4*ie+2]],
		       y_0[connectivity[4*ie+3]] };

      // create local element force accumulators
      float fxe[4] = { 0.0, 0.0, 0.0, 0.0 };
      float fye[4] = { 0.0, 0.0, 0.0, 0.0 };

      // loop over all quadrature points in the current element
      for (int q = 0; q < 4; q++) {
	// compute the shape function gradients w.r.t. parent coordinates
	float dSF_dxi[4]  = { 0.25*XIa[0]*(1.0+gauss*ETAa[0]*ETAa[q]),
			      0.25*XIa[1]*(1.0+gauss*ETAa[1]*ETAa[q]),
			      0.25*XIa[2]*(1.0+gauss*ETAa[2]*ETAa[q]),
			      0.25*XIa[3]*(1.0+gauss*ETAa[3]*ETAa[q]) };
	float dSF_deta[4] = { 0.25*ETAa[0]*(1.0+gauss*XIa[0]*XIa[q]),
			      0.25*ETAa[1]*(1.0+gauss*XIa[1]*XIa[q]),
			      0.25*ETAa[2]*(1.0+gauss*XIa[2]*XIa[q]),
			      0.25*ETAa[3]*(1.0+gauss*XIa[3]*XIa[q]) };

	// compute the Jacobian of the isoparametric transformation
	// to the element's reference coordinates
	float J[2][2] = { {x0e[0]*dSF_dxi[0]+x0e[1]*dSF_dxi[1]
			  +x0e[2]*dSF_dxi[2]+x0e[3]*dSF_dxi[3],
			   x0e[0]*dSF_deta[0]+x0e[1]*dSF_deta[1]
		          +x0e[2]*dSF_deta[2]+x0e[3]*dSF_deta[3]},
			  {y0e[0]*dSF_dxi[0]+y0e[1]*dSF_dxi[1]
	       	          +y0e[2]*dSF_dxi[2]+y0e[3]*dSF_dxi[3],
			   y0e[0]*dSF_deta[0]+y0e[1]*dSF_deta[1]
			  +y0e[2]*dSF_deta[2]+y0e[3]*dSF_deta[3]} };

	// compute the initial volume of the integration point
	// (the determinant of J), and its inverse
	float dV = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	float idetJ = 1.0 / dV;

	// compute the gradients w.r.t. the element's initial coordinates
	float grad_x0[4] = { (+dSF_dxi[0]*J[1][1]-dSF_deta[0]*J[1][0])*idetJ,
			     (+dSF_dxi[1]*J[1][1]-dSF_deta[1]*J[1][0])*idetJ,
			     (+dSF_dxi[2]*J[1][1]-dSF_deta[2]*J[1][0])*idetJ,
			     (+dSF_dxi[3]*J[1][1]-dSF_deta[3]*J[1][0])*idetJ };
	float grad_y0[4] = { (-dSF_dxi[0]*J[0][1]+dSF_deta[0]*J[0][0])*idetJ,
			     (-dSF_dxi[1]*J[0][1]+dSF_deta[1]*J[0][0])*idetJ,
			     (-dSF_dxi[2]*J[0][1]+dSF_deta[2]*J[0][0])*idetJ,
			     (-dSF_dxi[3]*J[0][1]+dSF_deta[3]*J[0][0])*idetJ };

	// compute the deformation gradient: dx/dX0
	float F[2][2] = { {xe[0]*grad_x0[0]+xe[1]*grad_x0[1]
			  +xe[2]*grad_x0[2]+xe[3]*grad_x0[3],
			   xe[0]*grad_y0[0]+xe[1]*grad_y0[1]
		          +xe[2]*grad_y0[2]+xe[3]*grad_y0[3]},
			  {ye[0]*grad_x0[0]+ye[1]*grad_x0[1]
			  +ye[2]*grad_x0[2]+ye[3]*grad_x0[3],
			   ye[0]*grad_y0[0]+ye[1]*grad_y0[1]
		          +ye[2]*grad_y0[2]+ye[3]*grad_y0[3]} };

	// compute the first Piola-Kirchhoff stress
	// (OpenACC currently does not support virtual methods!)
	model->update(F,&state_vars[n_state_vars*(4*ie+q)]);
	      
	// multiply by the differential volume dV
	float P_dV[2][2] = { {F[0][0]*dV, F[0][1]*dV},
			     {F[1][0]*dV, F[1][1]*dV} };

	// accumulate nodal forces
	fxe[0] += P_dV[0][0]*grad_x0[0] + P_dV[0][1]*grad_y0[0];
	fxe[1] += P_dV[0][0]*grad_x0[1] + P_dV[0][1]*grad_y0[1];
	fxe[2] += P_dV[0][0]*grad_x0[2] + P_dV[0][1]*grad_y0[2];
	fxe[3] += P_dV[0][0]*grad_x0[3] + P_dV[0][1]*grad_y0[3];
	fye[0] += P_dV[1][0]*grad_x0[0] + P_dV[1][1]*grad_y0[0];
	fye[1] += P_dV[1][0]*grad_x0[1] + P_dV[1][1]*grad_y0[1];
	fye[2] += P_dV[1][0]*grad_x0[2] + P_dV[1][1]*grad_y0[2];
	fye[3] += P_dV[1][0]*grad_x0[3] + P_dV[1][1]*grad_y0[3];
      } // for q = ...

      // store the average plot state
      plot[ie] = 0.0f;
      
      // scatter the accumulated element forces
      for (int a = 0; a < 4; a++) {
#pragma acc atomic update
	fx[connectivity[4*ie+a]] += fxe[a];
#pragma acc atomic update
	fy[connectivity[4*ie+a]] += fye[a];
      } // for a = ...
    } // for ie = ...

    // update nodal positions
#pragma acc parallel loop async(1) present(x,y,vx,vy,fx,fy,m) firstprivate(gx,gy,dt,idt)
    for (int in = 0; in < n_nodes; in++) {
      vx[in] += (-fx[in]*m[in] + gx)*dt;
      x[in]  += vx[in]*dt;
      vy[in] += (-fy[in]*m[in] + gy)*dt;
      y[in]  += vy[in]*dt;
      if (y[in] < 0.0) {
	float gap = y[in];
	y[in]  -= gap;
	vy[in] -= gap*idt;
      } // enforce bottom surface contact
      // re-initialize nodal forces
      fx[in] = 0.0;
      fy[in] = 0.0;
    } // for in = ...
    
  } // while (t < t_final)

#pragma acc wait(1)

  }

    auto t2 = std::chrono::high_resolution_clock::now();

  std::cout << "finished loop over time steps." << std::endl;

  std::cout << "simulation run-time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()
	    << " milliseconds" << std::endl;
  
} // TimeIntegrate()

void Update(void) {
  TimeIntegrate();
  glutPostRedisplay();
}

void InitGL(void) {
  glClearColor(0.9f,0.9f,0.9f,1);
  glEnable(GL_POINT_SMOOTH);
  glPointSize(H/2.f);
  glMatrixMode(GL_PROJECTION);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void Render(void) {
  glClear(GL_COLOR_BUFFER_BIT);
	
  glLoadIdentity();
  glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);

  glBegin(GL_QUADS);
  for (int ie = 0; ie < n_elems; ie++) {
    float blue = fmax(0.0,fmin(0.8,+plot[ie]/max_plot));
    float red  = fmax(0.0,fmin(0.8,-plot[ie]/max_plot));
    for (int a = 0; a < 4; a++) {
      glColor4f(0.8f-blue, 0.8f-red-blue, 0.8f-red, 1);
      glVertex2f(x[connectivity[4*ie+a]],y[connectivity[4*ie+a]]);
    } // for a = ...
  } // for ie = ...
  glEnd();

  glutSwapBuffers();
}

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y) {   
  switch(c) {
  case 'r': 
  case 'R':  
    InitFEM(); 
    break;
  }
}

void Arrows(int key, __attribute__((unused)) int x, __attribute__((unused)) int y) {   
  switch(key) {
  case GLUT_KEY_LEFT: 
    break;
  case GLUT_KEY_RIGHT:
    break;
  case GLUT_KEY_UP:
    break;
  case GLUT_KEY_DOWN:
    break;
  }
}

int main(int argc, char** argv) {
  glutInitWindowSize(WINDOW_WIDTH,WINDOW_HEIGHT);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutInit(&argc, argv);
  glutCreateWindow("FEM Library Demo");
  glutDisplayFunc(Render);
  glutIdleFunc(Update);
  glutKeyboardFunc(Keyboard);
  glutSpecialFunc(Arrows);

  InitGL();
  InitFEM();

  glutMainLoop();
  return 0;
} // main()
