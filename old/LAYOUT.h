#pragma once // compiler directive: only open this file once -- shortens compile times
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
class LAYOUT 
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    enum workaround1{d=TV::m};
public:
    T initial_time;
    int first_frame,last_frame;
    T frame_rate;
    int restart;
    std::string frame_title;
    int write_substeps_level;
    bool write_debug_data;
    std::string output_directory;

    T cfl;

    GRID<TV> mac_grid;
    THREAD_QUEUE* thread_queue;
    PROJECTION_COLLIDABLE_UNIFORM<GRID<TV> > projection;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<GRID<TV>,T, AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> > >,LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > > advection_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> boundary_scalar;
    BOUNDARY_UNIFORM<GRID<TV>,T> *boundary;
    T_LEVELSET levelset;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_boundary;
    RANGE<TV> source;
    pthread_mutex_t lock;

    const STREAM_TYPE stream_type;
//	const int number_of_frames; 	// Total number of frames
//    const T frame_time;          	// Frame (snapshot) interval
    T Time_At_Frame(const int frame) const;


    LAYOUT(const STREAM_TYPE stream_type_input);
	~LAYOUT();

    T CFL(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities);
    void CFL_Threaded(RANGE<TV_INT>& domain,ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,T& dt);

    void Initialize();
    T Maximum_Dt() const;
    void Initialize_Grid(TV_INT counts,RANGE<TV> domain)
        {mac_grid.Initialize(counts,domain,true);}
    virtual void Write_Output_Files(const int frame);
    virtual void Read_Output_Files(const int frame);
    virtual void Set_Boundary_Conditions(const T time);

    void Write_Output(const int frame) const;
};


}
