/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Heng Zhang, Ling Liu (Utah State University)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_potential.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "pair_lj_cut_coul_long.h"
#include "ewald.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

//test
#include <fstream>
#include <iostream>
#include <cstdio>
#include <iomanip>

using namespace std;
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

#define DANGER_ZONE    0.90
#define e  0.00000001

/* ---------------------------------------------------------------------- */

FixPot::FixPot(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{ 
  if (narg != 8) error->all(FLERR,"Illegal fix potential command");   

  // extract from pair
  cut_coul = cut_coulsq = 0.0;
  
  nevery     =  atoi(arg[3]);
  flag		 =  atoi(arg[4]);
  potentialp =  atof(arg[5])/2.;	
  potentialn = -atof(arg[5])/2.;
  
  tolerance  =  atof(arg[6]);
  pertype_parameters(arg[7]);	

  n = n_cap = N = nmax = m_fill = m_cap = 0;
  pack_flag = 0;		
  nprev = 0;
  s = b_s = ug = NULL;	
  unitk = NULL;
  kmax = kcount = 0;
  cs = sn = NULL;
  cstr = sstr = NULL;

  cut_coulsq = cut_coul*cut_coul;	
  A = AM = AMA = NULL;

  // register with Atom class
  grow_arrays(atom->nmax);			
  atom->add_callback(0);			

  for( int i = 0; i < atom->nmax; i++ )
  	 for (int j = 0; j < atom->nmax; ++j )
  		 A[i][j] = AM[i][j] = AMA[i][j] = 0;

  for (int i = 0; i < 2*atom->nmax ; i++)
	 for (int j = 0; j<3; j++ )
		for (int k = 0 ; k < atom->nmax ; k++ )
			cs[i][j][k] = sn[i][j][k] = 0; 

  for( int i = 0; i < 2*atom->nmax; i++ )
  	 for (int j = 0; j < atom->nmax; ++j )
  		 cstr[i][j] = sstr[i][j] = 0;
		 
}

/* ---------------------------------------------------------------------- */

FixPot::~FixPot()
{
  // unregister callbacks to this fix from Atom class
  atom->delete_callback(id,0);

  memory->destroy(A);
  memory->destroy(AM);
  memory->destroy(AMA);
  memory->destroy(cs);
  memory->destroy(sn);
  memory->destroy(cstr);
  memory->destroy(sstr);

  deallocate_storage();
  deallocate_matrix();
}

/* ---------------------------------------------------------------------- */

int FixPot::setmask() 
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPot::pertype_parameters(char *arg)
{
  if (!strcmp(arg,"lj/cut/coul/long")) {							
    Pair *pair = force->pair_match("lj/cut/coul/long",1);
    if (pair == NULL) error->all(FLERR,"No pair lj/cut/coul/long for fix potential");

    int tmp;
    double *cutoff = (double *) pair->extract("cut_coul",tmp);	

    if (cutoff == NULL)		
      error->all(FLERR,
                 "Fix potential could not cut-off from pair lj/cut/coul/long");

    cut_coul = *cutoff;		
    return;
  } else {
    error->all(FLERR,"No pair lj/cut/coul/long for fix potential");
  }
}

/* ---------------------------------------------------------------------- */

void FixPot::allocate_storage()
{
  nmax = atom->nmax;

  memory->create(s,nmax,"potential:s"); 
  memory->create(b_s,nmax,"potential:b_s");
  memory->create(ug,2*nmax,"potential:ug");
  memory->create(unitk,2*nmax,"potential:unitk");
}

/* ---------------------------------------------------------------------- */

void FixPot::deallocate_storage()
{
  memory->destroy( s );
  memory->destroy( b_s );
  memory->destroy( ug );
  memory->destroy( unitk );
}

/* ---------------------------------------------------------------------- */

void FixPot::reallocate_storage()
{
  deallocate_storage();
  allocate_storage();
  init_storage();
}

/* ---------------------------------------------------------------------- */

void FixPot::allocate_matrix()
{
  int i,ii;
  int mincap = 1000;
  double safezone = 10.0;

  n = atom->nlocal; 
  n_cap = MAX( (int)(n * safezone), mincap ); 

  int m = 0;								
  for( ii = 0; ii < list->inum; ii++ ) {	
    i = list->ilist[ii];					
    m += list->numneigh[i];					
  }
  m_cap = MAX( (int)(m * safezone), mincap);
}

/* ---------------------------------------------------------------------- */

void FixPot::deallocate_matrix()
{
}

/* ---------------------------------------------------------------------- */

void FixPot::reallocate_matrix()
{
  deallocate_matrix();
  allocate_matrix();
}

/* ---------------------------------------------------------------------- */

void FixPot::init()
{
  if (!atom->q_flag) error->all(FLERR,"Fix potential requires atom attribute q");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 1;	
  neighbor->requests[irequest]->full = 0;	
}

/* ---------------------------------------------------------------------- */

void FixPot::init_list(int id, NeighList *ptr) 
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixPot::setup_pre_force(int vflag)
{
  neighbor->build_one(list->index);	
									  
  deallocate_storage();				  
  allocate_storage();

  init_storage();

  deallocate_matrix();
  allocate_matrix();

  pre_force(vflag);			
}							
							
/* ---------------------------------------------------------------------- */

void FixPot::min_setup_pre_force(int vflag)
{
  setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPot::init_storage()			
{
  int itype;
  int *type = atom->type;
  double qqrd2e = force->qqrd2e;
  int k = 0;

  N = atom->nlocal + atom->nghost;	
  for( int i = 0; i < N; i++ ) {
	itype = type[i];
	s[i] = 0.0;
	b_s[i] = 0.0;
	if (itype > 4) {
	  if (itype == 5)	   b_s[k] = qqrd2e*potentialp;
	  else if (itype ==6 ) b_s[k] = qqrd2e*potentialn;
	  k++;
	}
  }
}

/* ---------------------------------------------------------------------- */

void FixPot::pre_force(int vflag)
{
  double t_start, t_end;
  int	 count;
  if (update->ntimestep % nevery) return;
  if ( comm->me == 0 ) t_start = MPI_Wtime();
  
  n = atom->nlocal;
  N = atom->nlocal + atom->nghost;
  
  int ghost = atom->nghost;
  if( atom->nmax > nmax ) reallocate_storage();
  if( n > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE ) 
    reallocate_matrix();

  /* ---------------------------------------------------------------------- */  

  xprd = domain->xprd;
  yprd = domain->yprd;
  zprd = domain->zprd;
  volume = xprd * yprd * zprd;
  
  // prepare unitk[]
  unitk[0] = 2.0*MY_PI/xprd;
  unitk[1] = 2.0*MY_PI/yprd;
  unitk[2] = 2.0*MY_PI/zprd;

  // prepare q2 and accuracy
  qsqsum = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    qsqsum += atom->q[i]*atom->q[i];
  }
  double tmp;
  MPI_Allreduce(&qsqsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsqsum = tmp;
  double q2 = qsqsum * force->qqrd2e / force->dielectric;

  double accuracy_relative = force->kspace->accuracy_relative;
  double two_charge_force  = force->kspace->two_charge_force;
  double accuracy = accuracy_relative * two_charge_force;

  // prepare kmax
  bigint natoms = atom->natoms;

  double err;
  kxmax = 1;
  kymax = 1;
  kzmax = 1;

  err = rms(kxmax,xprd,natoms,q2);
  while (err > accuracy) {
    kxmax++;
    err = rms(kxmax,xprd,natoms,q2);
  }

  err = rms(kymax,yprd,natoms,q2);
  while (err > accuracy) {
    kymax++;
    err = rms(kymax,yprd,natoms,q2);
  }

  err = rms(kzmax,zprd,natoms,q2);
  while (err > accuracy) {
    kzmax++;
    err = rms(kzmax,zprd,natoms,q2);
  }

  kmax = MAX(kxmax,kymax);
  kmax = MAX(kmax,kzmax);
  kmax3d = 4*kmax*kmax*kmax + 6*kmax*kmax + 3*kmax;

  // prepare gsqmx
  double gsqxmx = unitk[0]*unitk[0]*kxmax*kxmax;
  double gsqymx = unitk[1]*unitk[1]*kymax*kymax;
  double gsqzmx = unitk[2]*unitk[2]*kzmax*kzmax;
  gsqmx = MAX(gsqxmx,gsqymx);
  gsqmx = MAX(gsqmx,gsqzmx);

  coeffs();

  /* ---------------------------------------------------------------------- */
  
  count = init_matvec();
  //if (count != 0){	
  gaussian(count,AM,b_s,s);
  backup_Q(s);	
  //}
 
  force->kspace->setup();	
  if( comm->me == 0 ) {
	  t_end = MPI_Wtime();	
      potential_time = t_end - t_start;
  }
}

/* ---------------------------------------------------------------------- */

void FixPot::min_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

int FixPot::init_matvec()
{
  int ii, i, j, itype;
  int *type = atom->type;
  j = 0;

  for( ii = 0; ii < list->inum; ++ii ) {
    i = list->ilist[ii];
	itype = type[i];
	if (itype > 4) j++;
  }
  
  pack_flag = 1;
  comm->forward_comm_fix(this); //Dist_vector( s ); 
  // fill-in H matrix
  compute_ewald();
  
  return j;
}

/* ---------------------------------------------------------------------- */

void FixPot::compute_ewald() 
{
  int i, j, itype, jtype;
  double xtmp, ytmp, ztmp;
  double r, rsq, delr[3];
  int    *type = atom->type;
  double **x = atom->x;
  double *q = atom->q;
  double a = force->kspace->g_ewald;
  double erfcr;
  double qqrd2e = force->qqrd2e;
  int    buf[50000];
  double volume = domain->xprd*domain->yprd*domain->zprd;

  // initialization
  for (i = 0; i< atom->nlocal+atom->nghost ;i++ )
	for (j = 0; j < atom->nlocal+atom->nghost; j++ )
		AM[i][j] = 0.0;

  for (i=0;i<atom->nlocal;i++ ) {
	if (type[i]==5)			b_s[i] = qqrd2e*potentialp;
    else if (type[i]==6)	b_s[i] = qqrd2e*potentialn;
	else					b_s[i] = 0.0;
	if (flag == 1)			b_s[i] *= (rand()%100+1)*0.01;
  }

  // short-range term
  for ( i = 0 ; i < atom->nlocal+atom->nghost ; i++ ) {
	itype = type[i];
	if (itype < 5) continue;

	xtmp = x[i][0];
	ytmp = x[i][1];
	ztmp = x[i][2];

	for ( j = 0 ; j < atom->nlocal+atom->nghost ; j++ )
	{
		jtype = type[j];
		delr[0] = xtmp - x[j][0];
		delr[1] = ytmp - x[j][1];
		delr[2] = ztmp - x[j][2];
		rsq = delr[0]*delr[0] + delr[1]*delr[1] + delr[2]*delr[2];

		if (rsq > cut_coulsq) continue; 
		r = sqrt(rsq);
		erfcr = erfc(a*r);

		if (jtype > 4 && i != j)	AM[i][j]= 1./r*qqrd2e*erfcr;	
		else if (jtype < 5)			b_s[i] += -qqrd2e*q[j]/r*erfcr;
	}
  }
  
  printf("Short-range term complete.\n");

  // long-range term
  double self = compute_self();
  for ( i = 0 ; i < atom->nlocal ; i++ ) {
	if ( type[i] < 5 ) continue;
	// Ewald self term in H
	AM[i][i] = (-2*a/MY_PIS - MY_PI/(a*a*volume))*qqrd2e + self; 
	// non-diagonal term in H and b-term
	for ( j = 0 ; j < atom->nlocal ; j++ ) {
		if ( type[j] > 4 && i != j)
			AM[i][j] += -MY_PI/(a*a*volume)*qqrd2e + compute_cross(i,j);
		else if ( type[j] < 5)
			b_s[i] += q[j]*MY_PI/(a*a*volume)*qqrd2e + compute_b(i,j);	
	}
  }
	
  printf("Long-range term complete.\n"); 

  // map b_s matrix
  // comm->reverse_comm_fix( this );	

  // read mapping
  ifstream myfile;
  int list[50000], me[50000], met, mett, listt;
  int k = atom->nlocal;
  myfile.open("list1.txt");
  for (int i = 0 ; i < 20 /*35278 27697 14718 12821*/; i++ ) {
	  myfile >> met >> mett >> listt;
	  if (met == comm->me)
	  {	 
		  me[k] = met;
		  list[k] = listt; 
		  k++;
	  }
  }
  myfile.close();
    
  for (int i = atom->nlocal ; i < atom->nlocal+atom->nghost ; i++ ) {
	  while (list[i]>atom->nlocal-1) {
		  j = list[i];
		  list[i] = list[j];
	  } 
	  buf[i] = list[i];
  }

  for (int i=0 ; i < atom->nlocal ; i++ )
	for (int j=atom->nlocal; j < atom->nlocal+atom->nghost ; j++ )
	  AM[i][buf[j]] += AM[i][j];

  printf("mapping is done.\n");
  
  ofstream file;
  file.open("param.txt") ;
  file << atom->nlocal << endl;
  file << atom->nlocal+atom->nghost << endl;
  file << kcount << endl;
  file << force->kspace->g_ewald << endl;
  file << volume << endl;
  file.close();

  ofstream file2;
  file2.open("ug.txt");
  for (i=0; i < kcount; i++)
  	file2 << i << "\t" << ug[i] << endl;
  file2.close();

}

/* ---------------------------------------------------------------------- */

void FixPot::gaussian(int m, double **AB, double *bb, double *s)
{
	int itype, jtype, d1, d2, i, j;
	int *type = atom->type;
	double b[5000];
    	
	d1 = 1; 
	for ( int i = 0 ; i < n ; i++ ) {
		itype = type[i];
		if ( itype < 5 ) continue;
		d2 = 1;
		b[d1] = bb[i];
		for ( int j = 0 ; j < n ; j++) {
			jtype = type[j];
			if ( j == i ) {
				A[d1][d2] = AB[i][i];
				d2++;
			}
			if ( j != i && jtype > 4 ) {
				A[d1][d2] = AB[i][j];
				d2++;
			}
		}
		d1++;
	}

	ofstream file;
	file.open("A.txt");
	for (i=1; i<=m;i++ ) {
		for (j=1;j<=m ;j++ )
			file << A[i][j] << "\t";
		file << endl;
	}
	file.close();

	ofstream fil2;
	fil2.open("b.txt");
	for (i=1;i<=m;i++ )
		for (j=1;j<=m ;j++ )
			fil2 << b[i] << endl;
	fil2.close();

	// gauss elimination
	for(int k=1;k<m;k++){
        double ab_max=-1;
        int       max_ik;
        for(int i=k;i<=m;i++){
            if(abs(A[i][k])>ab_max){
                ab_max=abs(A[i][k]);
                max_ik=i;
            }
        }
        if(ab_max<e){
			printf("%f \t det A=0 \n",ab_max);
            break;
        }else if(max_ik!=k){
            double temp;
            for(int j=1;j<=m;j++){
                temp=A[max_ik][j];
                A[max_ik][j]=A[k][j];
                A[k][j]=temp;
            }
            temp=b[max_ik];
            b[max_ik]=b[k];
            b[k]=temp;
        }
        
        for(int i=k+1;i<=m;i++){
            A[i][k]/=A[k][k];
            for(int j=k+1;j<=m;j++){
                A[i][j]-=A[i][k]*A[k][j];
            }
            b[i]-=A[i][k]*b[k];
        }
        if(k<m-1)continue;
        else{
            if(abs(A[m][m])<e){
				printf("det A=0 \n");
                break;
            }else{
                s[m]=b[m]/A[m][m];
                for(int i=m-1;i>0;i--){
                    s[i]=b[i];
                    for(int j=i+1;j<=m;j++)
                        s[i]-=A[i][j]*s[j];
                    s[i]/=A[i][i];
                }
                
                //Print(m);
            }
        }
    }
   
}

/* ---------------------------------------------------------------------- */

void FixPot::Print(int m){
    printf("|-----------------------------\n");
    printf("|the results are:\n");
    for(int i=1;i<=m;i++){
        printf("s[%d]=  %lf\n",i,s[i]);
    }
    printf("|-----------------------------\n");
}

/* ---------------------------------------------------------------------- */

void FixPot::backup_Q(double *s)		
{
  int i, k, itype;
  int *type = atom->type;
  double *q = atom->q;
  double **x = atom->x;
  double delx,dely,delz;
    
  // update q and backup
  k = 1; 
  for( i = 0; i < n; ++i ) {
	  itype = type[i];
	  if (itype > 4) {
		  q[i] = s[k]; 
		  k++;
	  }
  }

  pack_flag = 2;
  comm->forward_comm_fix( this ); //Dist_vector( atom->q );
  
  ofstream file;
  file.open("atoms.txt");
  for(i=0; i< atom->nlocal+atom->nghost;i++)
  	file << i << "\t" << atom->x[i][0] << "\t" << atom->x[i][1] << "\t"
  				<< atom->x[i][2] << "\t" << q[i] << endl;
  file.close();
  	
}

/* ---------------------------------------------------------------------- */

int FixPot::pack_comm(int n, int *list, double *buf,
                          int pbc_flag, int *pbc)
{
  int m;
  ofstream myfile;
  myfile.open("list1.txt",std::fstream::app);
  if( pack_flag == 1) {
    for(m = 0; m < n; m++) { 
		buf[m] = s[list[m]];
		myfile << comm->me << "\t" << m << "\t" << list[m] << endl;
	}
  }
  else if( pack_flag == 2 )
    for(m = 0; m < n; m++) buf[m] = atom->q[list[m]];
  return 1;
  myfile.close();
}

/* ---------------------------------------------------------------------- */

void FixPot::unpack_comm(int n, int first, double *buf)
{
  int i, m;
  ofstream myfile;
  myfile.open("list2.txt",std::fstream::app);
  if( pack_flag == 1) {
    for(m = 0, i = first; m < n; m++, i++) {
		s[i] = buf[m];
		myfile << comm->me << "\t" << i << endl;
	}
  }
  else if( pack_flag == 2)
    for(m = 0, i = first; m < n; m++, i++) atom->q[i] = buf[m];
	myfile.close();
}

/* ---------------------------------------------------------------------- */

int FixPot::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m;
  for(m = 0, i = first; m < n; m++, i++) buf[m] = b_s[i]; 
  return 1; 
}

/* ---------------------------------------------------------------------- */

void FixPot::unpack_reverse_comm(int n, int *list, double *buf)
{
  for(int m = 0; m < n; m++) b_s[list[m]] += buf[m];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPot::memory_usage()
{
  double bytes;

  bytes += atom->nmax*11 * sizeof(double); // storage
  bytes += n_cap*2 * sizeof(int); // matrix...
  bytes += m_cap * sizeof(int);
  bytes += m_cap * sizeof(double);

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate H and arrays
------------------------------------------------------------------------- */

void FixPot::grow_arrays(int nmax)
{
  memory->grow(A,nmax,nmax,"potential:A");
  memory->grow(AM,nmax,nmax,"potential:AM");
  memory->grow(AMA,nmax,nmax,"potential:AMA");
  memory->grow(ug,2*nmax,"potential:ug");
  memory->grow(unitk,2*nmax,"potential:unitk");
  memory->grow(cs,2*nmax,3,nmax,"potential:cs");
  memory->grow(sn,2*nmax,3,nmax,"potential:sn");
  memory->grow(cstr,2*nmax,nmax,"potential:cstr");
  memory->grow(sstr,2*nmax,nmax,"potential:sstr");
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
-------------------------------------------------------------------------*/

int FixPot::pack_exchange(int i, double *buf) 
{ 
  return nprev;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPot::unpack_exchange(int nlocal, double *buf) 
{
  return nprev;
}

/* ---------------------------------------------------------------------- */

void FixPot::coeffs()
{
  int i,j,k,l,m,n,ic;
  double sqk,clpm,slpm;

  double **x = atom->x;
  int nlocal = atom->nlocal;

  n = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (ic = 0; ic < 3; ic++) {
    sqk = unitk[ic]*unitk[ic];
	if (sqk <= gsqmx) {
      for (i = 0; i < nlocal; i++) {
        cs[0][ic][i] = 1.0;
        sn[0][ic][i] = 0.0;
        cs[1][ic][i] = cos(unitk[ic]*x[i][ic]);
        sn[1][ic][i] = sin(unitk[ic]*x[i][ic]);
        cstr[n][i] = cs[1][ic][i];
        sstr[n][i] = sn[1][ic][i];
      }
	  n++;
    }
  }

  for (m = 2; m <= kmax; m++) {
    for (ic = 0; ic < 3; ic++) {
      sqk = m*unitk[ic] * m*unitk[ic];
      if (sqk <= gsqmx) {
        for (i = 0; i < nlocal; i++) {
          cs[m][ic][i] = cs[m-1][ic][i]*cs[1][ic][i] -
            sn[m-1][ic][i]*sn[1][ic][i];
          sn[m][ic][i] = sn[m-1][ic][i]*cs[1][ic][i] +
            cs[m-1][ic][i]*sn[1][ic][i];
          cstr[n][i] = cs[m][ic][i]; 
          sstr[n][i] = sn[m][ic][i]; 
        }
        n++;
      }
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]);
      if (sqk <= gsqmx) {
        for (i = 0; i < nlocal; i++) {
          cstr[n][i] = cs[k][0][i]*cs[l][1][i] - sn[k][0][i]*sn[l][1][i];
          sstr[n][i] = sn[k][0][i]*cs[l][1][i] + cs[k][0][i]*sn[l][1][i];
		}
		n++;
		for (i = 0; i < nlocal; i++) {
          cstr[n][i] = cs[k][0][i]*cs[l][1][i] + sn[k][0][i]*sn[l][1][i];
          sstr[n][i] = sn[k][0][i]*cs[l][1][i] - cs[k][0][i]*sn[l][1][i];
        }
		n++;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (l*unitk[1] * l*unitk[1]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
        for (i = 0; i < nlocal; i++) {
          cstr[n][i] = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
          sstr[n][i] = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
		}
		n++;
		for (i = 0; i < nlocal; i++) {
          cstr[n][i] = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
          sstr[n][i] = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
        }
        n++;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
        for (i = 0; i < nlocal; i++) {
          cstr[n][i] = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
          sstr[n][i] = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
		}
		n++;
		for (i = 0; i < nlocal; i++) {
          cstr[n][i] = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
          sstr[n][i] = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
        }
        n++;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]) +
          (m*unitk[2] * m*unitk[2]);
        if (sqk <= gsqmx) {
          for (i = 0; i < nlocal; i++) {
            clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
            slpm = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
            cstr[n][i] = cs[k][0][i]*clpm - sn[k][0][i]*slpm;
            sstr[n][i] = sn[k][0][i]*clpm + cs[k][0][i]*slpm;
		  }
		  n++;
		  for (i = 0; i < nlocal; i++) {
            clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
            slpm = -sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
            cstr[n][i] = cs[k][0][i]*clpm - sn[k][0][i]*slpm;
            sstr[n][i] = sn[k][0][i]*clpm + cs[k][0][i]*slpm;
		  }
		  n++;
		  for (i = 0; i < nlocal; i++) {
            clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
            slpm = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
            cstr[n][i] = cs[k][0][i]*clpm - sn[k][0][i]*slpm;
            sstr[n][i] = sn[k][0][i]*clpm + cs[k][0][i]*slpm;
		  }
		  n++;
		  for (i = 0; i < nlocal; i++) {
            clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
            slpm = -sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
            cstr[n][i] = cs[k][0][i]*clpm - sn[k][0][i]*slpm;
            sstr[n][i] = sn[k][0][i]*clpm + cs[k][0][i]*slpm;
          }
		  n++;
        }
      }
    }
  }

  ofstream file;
  file.open("cossin.txt");
  for (i=0; i< n ; i++)
  	for (j=0; j< atom->nlocal; j++)
  		file << i << "\t" << cstr[i][j] << "\t" << sstr[i][j] << endl;	
  file.close();
}

/* ----------------------------------------------------------------------
   compute self term long-range
------------------------------------------------------------------------- */

double FixPot::compute_self()
{
  int i,k,l,m;
  double sqk, self;
  double g_ewald = force->kspace->g_ewald;
  double g_ewald_sq_inv = 1.0 / (g_ewald*g_ewald);
  double preu = 4.0*MY_PI/volume;
  double qqrd2e = force->qqrd2e; 

  kcount = 0;
  // (k,0,0), (0,l,0), (0,0,m)

  for (m = 1; m <= kmax; m++) {
    sqk = (m*unitk[0]) * (m*unitk[0]);
    if (sqk <= gsqmx) {
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      kcount++;
    }
    sqk = (m*unitk[1]) * (m*unitk[1]);
    if (sqk <= gsqmx) {
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      kcount++;
    }
    sqk = (m*unitk[2]) * (m*unitk[2]);
    if (sqk <= gsqmx) {
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      kcount++;
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      sqk = (unitk[0]*k) * (unitk[0]*k) + (unitk[1]*l) * (unitk[1]*l);
      if (sqk <= gsqmx) {
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;

        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kymax; l++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[1]*l) * (unitk[1]*l) + (unitk[2]*m) * (unitk[2]*m);
      if (sqk <= gsqmx) {
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;

        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kxmax; k++) {
    for (m = 1; m <= kzmax; m++) {
      sqk = (unitk[0]*k) * (unitk[0]*k) + (unitk[2]*m) * (unitk[2]*m);
      if (sqk <= gsqmx) {
        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;

        ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
        kcount++;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kxmax; k++) {
    for (l = 1; l <= kymax; l++) {
      for (m = 1; m <= kzmax; m++) {
        sqk = (unitk[0]*k) * (unitk[0]*k) + (unitk[1]*l) * (unitk[1]*l) +
          (unitk[2]*m) * (unitk[2]*m);
        if (sqk <= gsqmx) {
          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;

          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;

          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;

          ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
          kcount++;
        }
      }
    }
  }

  self = 0.0;
  for (k = 0; k < kcount ; k++ ) self += ug[k];
  return self*qqrd2e;
  
}

/* ----------------------------------------------------------------------
   compute RMS accuracy for a dimension
------------------------------------------------------------------------- */ 

double FixPot::rms(int km, double prd, bigint natoms, double q2)
{
  double g_ewald = force->kspace->g_ewald; 

  double value = 2.0*q2*g_ewald/prd *
    sqrt(1.0/(MY_PI*km*natoms)) *
    exp(-MY_PI*MY_PI*km*km/(g_ewald*g_ewald*prd*prd));

  return value;
}

/* ---------------------------------------------------------------------- */

double FixPot::compute_cross(int i, int j)
{
  double cross = 0.0;
  double qqrd2e = force->qqrd2e;

  for (int k = 0; k < kcount ; k++ ) 
	cross += ug[k]*(cstr[k][i]*cstr[k][j] + sstr[k][i]*sstr[k][j]);  

  return cross*qqrd2e;
}

/* ---------------------------------------------------------------------- */

double FixPot::compute_b(int i, int j)
{
  double b_recip = 0.0;
  double qqrd2e = force->qqrd2e;

  for (int k = 0; k < kcount ; k++ ) 
	 b_recip += -ug[k]*atom->q[j]*
	  (cstr[k][i]*cstr[k][j] + sstr[k][i]*sstr[k][j]); 

  return b_recip*qqrd2e;
}