//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_THREADED.h>
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"

using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer, const std::string& title, int substep, int level)
{
    ((WATER_DRIVER<TV>*)writer)->Write_Substep(title, substep, level);
}
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
WATER_DRIVER(WATER_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this, &Write_Substep_Helper<TV>);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> WATER_DRIVER<TV>::
~WATER_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Execute_Main_Program
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if (example.restart) {current_frame=example.restart;}else {current_frame=example.first_frame;}    
    time=example.Time_At_Frame(current_frame);

    // mpi
    if (example.mpi_grid) 
		{example.mpi_grid->Initialize(example.domain_boundary);}
    example.projection.elliptic_solver->mpi_grid = example.mpi_grid;
    if(example.mpi_grid) 
		{example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);}
    else if (example.thread_queue) 
		{example.boundary=new BOUNDARY_THREADED<GRID<TV> >(*example.thread_queue,example.boundary_scalar);}    
    else 
		{example.boundary=&example.boundary_scalar;}

    //threading
    example.projection.elliptic_solver->thread_queue=example.thread_queue;

    // setup grids and velocities
    example.projection.Initialize_Grid(example.mac_grid);
    example.face_velocities.Resize(example.mac_grid);
    example.levelset.phi.Resize(example.mac_grid.Domain_Indices(3));
	example.temperature.Resize(example.mac_grid.Domain_Indices(3));		// Mike 05/03/14
    example.Initialize_Fields();
	

    // setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance(1e-9);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    example.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations=40;
    example.projection.collidable_solver->Use_External_Level_Set(example.levelset);    

    if(example.restart) example.Read_Output_Files(example.restart);
    
    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;
    constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.Set_Boundary_Conditions(time); // get so CFL is correct

    if(!example.restart) 
		{Write_Output_Files(example.first_frame);}
    output_number=example.first_frame;
}
//#####################################################################
// Scalar_Advance
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Scalar_Advance(const T dt,const T time)
{
    example.Get_Scalar_Field_Sources(time);
	
	// advect the levelset
    ARRAY<T,TV_INT> phi_ghost(example.mac_grid.Domain_Indices(3));
    example.boundary->Set_Fixed_Boundary(true,0);
    example.boundary->Fill_Ghost_Cells(example.mac_grid,example.levelset.phi,phi_ghost,dt,time,3);
    example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid,example.levelset.phi,phi_ghost,example.face_velocities,*example.boundary,dt,time); 
    example.boundary->Set_Fixed_Boundary(false);
    example.levelset.Fast_Marching_Method();

	// advect the temperatures -- Mike 5/3/2014
    ARRAY<T,TV_INT> temperature_ghost(example.mac_grid.Domain_Indices(3));
    example.boundary->Fill_Ghost_Cells(example.mac_grid,example.temperature,temperature_ghost,dt,time,3);
    example.advection_scalar.Update_Advection_Equation_Cell(example.mac_grid,example.temperature,temperature_ghost,example.face_velocities,*example.boundary,dt,time); 
}
//#####################################################################
// Convect
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Convect(const T dt,const T time)
{	
	// advect the velocities
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(example.mac_grid,3,false);
    example.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,time,3);
    example.advection_scalar.Update_Advection_Equation_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,face_velocities_ghost,*example.boundary,dt,time);
}
//#####################################################################
// Add_Body_Forces on may 2 2014 -- updated by mike on may 3 2014
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Add_Body_Forces(const T dt, const T time)
{
	T density1 = 1.0;
	T density2 = 1.1;

	// incorporate buoyancy
	// iterate through all cells and adjust velocity based on cell temperature, which liquid, and density
	
	


	T gravity = 0;// /(10*time+1);
	T buoyancy;

	//This is where we iterate through all 'water' cells and adjust the temperature
	// 'air' cells will be set to a constant value. 

    //For all dimensions: x, y, and maybe z
    for(int axis=1;axis<=TV::dimension;axis++)
    {
		// Look at both sides of this cell along this dimension (e.g. Left/Right, Up/Down, Front/Back)
		for(int axis_side=1;axis_side<=TV::dimension;axis_side++)
		{
			//Give this side a unique number relative to this cell (1-4 for 2d) or (1-6 for 3d)
			int side=2*(axis-1)+axis_side; //making a value for side...
			
			//Get the offset to get the index of the cell relative to the index of this face
			TV_INT interior_cell_offset;
			if(axis_side==1)
			{
				interior_cell_offset = TV_INT();
			}
			else
			{
				interior_cell_offset = -TV_INT::Axis_Vector(axis);
			} 
	
			//Step through every face
			for(typename GRID<TV>::FACE_ITERATOR iterator(example.mac_grid); iterator.Valid(); iterator.Next())
			{
					//Grab the index using the offset we calculated above
			        TV_INT cell=iterator.Face_Index()+interior_cell_offset;
				//std::cout<<"this is axis side = "<<axis_side;
				//std::cout<<", and this is cell = "<<cell<<std::endl;
				//std::cout<<", and this is cell x?= "<<cell.x<<std::endl; this returns the x component of this vector, and same for y below. -DTR 05/11/2014
				//std::cout<<", and this is cell y?= "<<cell.y<<std::endl;
			        if(example.levelset.phi(cell)<=0 && (axis ==1 || axis == 3))//If this cell is in water
				{
					//Calculate Buoyancy based on temperature
					buoyancy = ((example.resolution)/2.0 - cell.x)/(time*1+1);
				}
			        else //In air
				{
					//Set Buoyancy to a constant, since we are in the 'air'
					buoyancy = 0;
				}

					//We now have gravity & buoyancy. We will combine them and then set velocity based on all external forces. 
				T externalForces = gravity + buoyancy;

				//int axis=iterator.Axis();
				if (axis == 1) 
				{
		            example.face_velocities.Component(axis)(iterator.Face_Index()) += dt*externalForces; //sneaky gravity was 8.8 put to zero for test
			//		continue;
				}
				else if (axis == 3) 
				{
		            example.face_velocities.Component(axis)(iterator.Face_Index()) += (dt*externalForces); //sneaky gravity was 8.8 put to zero for test
			//		continue;
				}
		            //example.face_velocities.Component(axis)(iterator.Face_Index()) -= dt*externalForces; //sneaky gravity was 9.8 put to zero for test
			}
		}
	}
/* 
	for (typename GRID<TV>::CELL_ITERATOR i(mac_grid); i.Valid(); i.Next()) {
		product = i.Location()-center;
		product *= product;
		d = negative_radius_squared;
		for(int i=1; i<=TV::dimension; i++) d+=product(i);
		levelset.phi(i.Cell_Index()) = d;
	}*/


	




}
//#####################################################################
// Project
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Project(const T dt,const T time)
{
    example.Set_Boundary_Conditions(time+dt);
    example.projection.p*=dt; // rescale pressure for guess
    example.boundary->Apply_Boundary_Condition_Face(example.mac_grid,example.face_velocities,time+dt);
    example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method();
    example.projection.Make_Divergence_Free(example.face_velocities,dt,time);
    example.projection.collidable_solver->Set_Up_Second_Order_Cut_Cell_Method(false);
    example.projection.p*=(1/dt); // unscale pressure

//-DTR this code between here and the below comment can be commented out. we are not sure exactly what it does, however, it gives results either way.
    const int ghost_cells=7;    
    T delta=20*example.mac_grid.dX.Max();
    ARRAY<T,TV_INT> phi_ghost(example.mac_grid.Domain_Indices(3));
    example.boundary->Fill_Ghost_Cells(example.mac_grid,example.levelset.phi,phi_ghost,dt,time,3);
    for(int axis=1;axis<=GRID<TV>::dimension;axis++)
    {
        GRID<TV> face_grid=example.mac_grid.Get_Face_Grid(axis);
	ARRAY<T,TV_INT> phi_face(face_grid.Domain_Indices(),false);
	T_ARRAYS_BASE& face_velocity=example.face_velocities.Component(axis);
        ARRAY<bool,TV_INT> fixed_face(face_grid.Domain_Indices());
        for(typename GRID<TV>::FACE_ITERATOR iterator(example.mac_grid,0,GRID<TV>::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next())
	{
            TV_INT index=iterator.Face_Index();
	    phi_face(index)=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
            if(phi_face(index)<=0) fixed_face(index)=true;
	    if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0.0;//where phi > delta, where delta is 3*dx, as noted above, set these velocities = 0; 
	}
        LOG::cout<<"something..."<<std::endl;  // TODO(jontg): If this log statement doesn't appear, the code crashes in release mode...
        T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,ghost_cells);
	extrapolate.Set_Band_Width(3);
	extrapolate.Set_Custom_Seed_Done(&fixed_face);
        extrapolate.Extrapolate();
    }
 // this whole section of code was commented out, seemingly this has not caused an error and allows us to say that these two liquids are liquid and not just liquid and air. -DTR 05/11/2014
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
	for (int substep=1; !done; substep++) {
        // choose time step
        LOG::Time("Substep");
        T dt = example.cfl*example.CFL(example.face_velocities);        
        if (example.mpi_grid) 
			{example.mpi_grid->Synchronize_Dt(dt);}
        if (time+dt >= target_time) 
			{dt = target_time-time; done=true;}
        else if (time+2*dt >= target_time)
			{dt = .5*(target_time-time);}

	// advance the simulation
        Scalar_Advance(dt,time);
        Convect(dt,time);
        Add_Body_Forces(dt,time);
        Project(dt, time);
        time += dt;
	}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        LOG::cout<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame == example.first_frame) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
template class WATER_DRIVER<VECTOR<float,1> >;
template class WATER_DRIVER<VECTOR<float,2> >;
template class WATER_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class WATER_DRIVER<VECTOR<double,1> >;
template class WATER_DRIVER<VECTOR<double,2> >;
template class WATER_DRIVER<VECTOR<double,3> >;
#endif
