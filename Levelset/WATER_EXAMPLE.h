//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WATER_EXAMPLE__
#define __WATER_EXAMPLE__
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/PROJECTION_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
namespace PhysBAM{

template<class TV>
class WATER_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    enum workaround1{d=TV::m};

public:
    STREAM_TYPE stream_type;
    T initial_time;
    int first_frame;
	int last_frame;
    T frame_rate;
    int restart;
    std::string frame_title;
    int write_substeps_level;
    bool write_debug_data;
    std::string output_directory;
    int height;  //hey mike, we added this-DA 05/02/2014
    int resolution;  //hey mike, we added this-DA 05/02/2014
    T cfl;

    GRID<TV> mac_grid;
    MPI_UNIFORM_GRID<GRID<TV> > *mpi_grid;
    THREAD_QUEUE* thread_queue;    
    PROJECTION_COLLIDABLE_UNIFORM<GRID<TV> > projection;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
	
    //PROJECTION_COLLIDABLE_UNIFORM<GRID<TV> > cell_heat;// we are assuming this makes an array of type double which is length of cell_index long, we think it has one entry per grid cell.
    //double cell_heat[FACE_INDEX<TV::dimension>];
	
	// I think all we need to keep track of the heat / temperature is just an arry -- Mike 5/3/2014		
	ARRAY<T,TV_INT> temperature;

    ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>,T, AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> >  >,LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > > advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary;
    T_LEVELSET levelset;
    VECTOR<VECTOR<bool, 2>, TV::dimension> domain_boundary;    
    //RANGE<TV> source;
    pthread_mutex_t lock;

    WATER_EXAMPLE(const STREAM_TYPE stream_type_input,int refine=0);
    virtual ~WATER_EXAMPLE();

    T CFL(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities);
    void CFL_Threaded(RANGE<TV_INT>& domain, ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities, T& dt);
    T Time_At_Frame(const int frame) const {return initial_time+(frame-first_frame)/frame_rate;}
	void Initialize_Grid(TV_INT counts, RANGE<TV> domain) {mac_grid.Initialize(counts, domain, true);}
	void Initialize_Fields(); //initializes level set phi and face velocities
	void Get_Scalar_Field_Sources(const T time_velocities);

    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Set_Boundary_Conditions(const T time);

//#####################################################################
};
}
#endif
