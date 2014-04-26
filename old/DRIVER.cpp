#include "DRIVER.h"
#include "LAYOUT.h"

using namespace PhysBAM;

// constructor
template<class TV>
DRIVER<TV>::DRIVER(LAYOUT<TV>& layout) 
	:layout(layout) {}

// destructor
template<class TV>
DRIVER<TV>::~DRIVER() {}

// Run
template<class TV>
void DRIVER<TV>::Run()
{
	Initialize();
	layout.Write_Output(0);// this outputs the grid, -DA 04/26/2014
	time=0;

	for(int frame=1; frame<=layout.last_frame; frame++) {//this simulates the frames and outputs them -DA 04262014
		Simulate_Frame(frame);
		layout.Write_Output(frame);
	}
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Initialize()
{
    //DEBUG_SUBSTEPS::Set_Write_Substeps_Level(layout.write_substeps_level);

    // setup time
    if(layout.restart) current_frame=layout.restart;
    else current_frame=layout.first_frame;
    
    time=layout.Time_At_Frame(current_frame);

    // set layout.boundary
    layout.boundary=&layout.boundary_scalar;

    // setup grids and velocities
    layout.projection.Initialize_Grid(layout.mac_grid);
    layout.face_velocities.Resize(layout.mac_grid);
    layout.levelset.phi.Resize(layout.mac_grid.Domain_Indices(3));
    layout.Initialize_Fields();

    // setup laplace
    layout.projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
    layout.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    layout.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    layout.projection.elliptic_solver->pcg.cg_restart_iterations=40;
    layout.projection.collidable_solver->Use_External_Level_Set(layout.levelset);

    if(layout.restart) layout.Read_Output_Files(layout.restart);

    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    layout.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    layout.Set_Boundary_Conditions(time); // get so CFL is correct

    if(!layout.restart) Write_Output_Files(layout.first_frame);
    output_number=layout.first_frame;
}
//#####################################################################

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


