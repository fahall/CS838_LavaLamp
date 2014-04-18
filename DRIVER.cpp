#include "DRIVER.h"
#include "LAYOUT.h"

using namespace PhysBAM;

// constructor
template<class TV>
DRIVER<TV>::DRIVER(LAYOUT<TV>& layout_input) 
	:layout(layout_input) {}

// destructor
template<class TV>
DRIVER<TV>::~DRIVER() {}

// Run
template<class TV>
void DRIVER<TV>::Run()
{
	layout.Initialize();
	layout.Write_Output(0);
	time=0;

	for(int frame=1; frame<=layout.number_of_frames; frame++) {
		Simulate_Frame(frame);
		layout.Write_Output(frame);
	}
}

// Simulate_Frame : TAKE TIME STEPS TO ADVANCE TO THE NEXT FRAME
template<class TV>
void DRIVER<TV>::Simulate_Frame(const int frame)
{
    T frame_end_time = layout.frame_time * (T)frame;
	T dt_max;
	T dt;

    for(;time<frame_end_time;time+=dt){
        dt_max = layout.Maximum_Dt();
        dt = std::min(dt_max, (T)1.001*(frame_end_time - time));
        Simulate_Time_Step(time,dt);
	}
}


// Simulate_Time_Step : ADVANCE ONE STEP FORWARD IN TIME
template<class TV>
void DRIVER<TV>::Simulate_Time_Step(const T time,const T dt)
{
 // #### TODO ####

}

template class DRIVER<VECTOR<float,2> >;
template class DRIVER<VECTOR<float,3> >;


