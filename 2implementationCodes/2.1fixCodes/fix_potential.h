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

#ifdef FIX_CLASS

FixStyle(potential,FixPot)

#else

#ifndef LMP_FIX_POTENTIAL_H
#define LMP_FIX_POTENTIAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPot : public Fix {
 public:
  FixPot(class LAMMPS *, int, char **);
  ~FixPot();
  int setmask();
  void init();
  void init_list(int,class NeighList *);
  void init_storage();
  void setup_pre_force(int);
  void pre_force(int);

  void min_setup_pre_force(int);
  void min_pre_force(int);

  double potential_time;

 private:
  int nevery,pairflag,flag;
  int n, N, m_fill;
  int n_cap, nmax, m_cap;
  int pack_flag;
  class NeighList *list;
  int mask;

  double tolerance;     // tolerance for the norm of the rel residual in CG
  double cut_coul,cut_coulsq;	

  double *s, *b_s,*ug; 
  double **A, **AM, **AMA;
  int	 nprev;
  double potentialp, potentialn;
  double *x;

  // long-range term
  double *unitk;
  double qsqsum, gsqmx, alpha;
  int kxmax, kymax, kzmax, kmax, kmax3d;
  int kcount;
  double xprd, yprd, zprd, volume; 
  double ***cs,***sn;
  double **cstr, **sstr;

  void pertype_parameters(char*);
  void allocate_storage();
  void deallocate_storage();
  void reallocate_storage();
  void allocate_matrix();
  void deallocate_matrix();
  void reallocate_matrix();

  double rms(int, double, bigint, double);
  void	 coeffs();
  void   compute_ewald();
  int    init_matvec();
  double compute_self();
  double compute_cross(int, int);
  double compute_b(int, int);
  void gaussian(int, double **, double *,double *);
  void backup_Q(double *);

  void Print(int);
  //void PrintProc(int, double *);
  int  pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int  pack_reverse_comm(int, int, double *);     
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();
  void grow_arrays(int);
  int  pack_exchange(int, double *);
  int  unpack_exchange(int, double *);
  
  class Ewald *ewald;

};

}

#endif
#endif
