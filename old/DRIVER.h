
#pragma once // compiler directive: only open this file once -- shortens compile times
/*
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>

namespace PhysBAM{

// I think this line makes sure the LAYOUT class is compiled first
template<class TV> class LAYOUT;

template<class TV>
class DRIVER
{
	typedef typename TV::SCALAR T;

	public:
	LAYOUT<TV>& layout;
	T time;
    
	DRIVER(LAYOUT<TV>& layout_input);
	~DRIVER();
	void Run();
	void Simulate_Frame(const int frame);
	void Simulate_Time_Step(const T time,const T dt);
};
}*/
#ifndef __WATER_DRIVER__
#define __WATER_DRIVER__
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
namespace PhysBAM{


template<class TV> class LAYOUT;

template<class TV>
class DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef EXTRAPOLATION_UNIFORM<GRID<TV>,T> T_EXTRAPOLATION_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_BASE T_ARRAYS_BASE;

protected:
    int current_frame;
    T time;
    int output_number;

    LAYOUT<TV>& example;
public:
    DRIVER(LAYOUT<TV>& layout);
    virtual ~DRIVER();

    void Scalar_Advance(const T dt,const T time);
    void Convect(const T dt,const T time);
    void Add_Forces(const T dt,const T time);
    void Project(const T dt,const T time);
    void Initialize();
    void Advance_To_Target_Time(const T target_time);
    void Simulate_To_Frame(const int frame_input);
    void Write_Output_Files(const int frame);
    void Write_Substep(const std::string& title,const int substep,const int level=0);

//#####################################################################
};
}
#endif
