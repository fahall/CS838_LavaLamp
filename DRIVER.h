
#pragma once // compiler directive: only open this file once -- shortens compile times

#include <PhysBAM_Tools/Vectors/VECTOR.h>

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
}
