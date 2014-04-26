#include "LAYOUT.h"
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>

using namespace PhysBAM;
//#####################################################################
// LAYOUT Constructor
//#####################################################################
template<class TV>
LAYOUT<TV>::LAYOUT(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    restart(0),write_debug_data(false),output_directory("output"),cfl(.9),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    projection(mac_grid,false,false,thread_queue),advection_scalar(thread_queue),
levelset(mac_grid,*new ARRAY<T,TV_INT>())
{}
//#####################################################################
// CFL or courant number, determines size of time step due to grid spacing
//#####################################################################
template<class TV> typename TV::SCALAR LAYOUT<TV>::
CFL(ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities)
{
    T dt=FLT_MAX;
    DOMAIN_ITERATOR_THREADED_ALPHA<LAYOUT<TV>,TV>(mac_grid.Domain_Indices(),thread_queue).template Run<ARRAY<T,FACE_INDEX<TV::dimension> >&,T&>(*this,&LAYOUT::CFL_Threaded,face_velocities,dt);
    return dt;
}
//#####################################################################
// LAYOUT destructor, empty because stub function, not needed
//#####################################################################
template<class TV>
LAYOUT<TV>::~LAYOUT()
{}
//#####################################################################
// get time at input frame based on frame rate
//#####################################################################
T Time_At_Frame(const int frame) const
{return initial_time+(frame-first_frame)/frame_rate;}

template<class TV>
void LAYOUT<TV>::Initialize()
{}


template<class TV>
typename TV::SCALAR LAYOUT<TV>::Maximum_Dt() const
{
    return 1e-3;
}
//#####################################################################
// Set_Boundary_Conditions
//#####################################################################
template<class TV> void LAYOUT<TV>::
Set_Boundary_Conditions(const T time)
{
    projection.elliptic_solver->psi_D.Fill(false);
    projection.elliptic_solver->psi_N.Fill(false);

    for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++)
    {
        int side=2*(axis-1)+axis_side;//making a value for side...
//if axis_side==1 then interior_cell_offset=TV_INT() else -TV_INT
        TV_INT interior_cell_offset, exterior_cell_offset, boundary_face_offset;
        if(axis_side==1)
        {
                interior_cell_offset = TV_INT();
                exterior_cell_offset = -TV_INT::Axis_Vector(axis);
                boundary_face_offset = TV_INT::Axis_Vector(axis);
        }
        else
        {
                interior_cell_offset = -TV_INT::Axis_Vector(axis);
                exterior_cell_offset = TV_INT();
                boundary_face_offset = -TV_INT::Axis_Vector(axis);
        }

        if(domain_boundary(axis)(axis_side))
        {//this is a boundary
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next())
            {
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                if(levelset.phi(face+interior_cell_offset)<=0)//checks if level set phi is less than or equal to zero(not water) at boundary face
                {

                    if(face_velocities.Component(axis).Valid_Index(face))//is this a valid position for face velocities? if so
                        {
                        projection.elliptic_solver->psi_N.Component(axis)(face)=true;//pressure solver is looking for a valid normal component
                        face_velocities.Component(axis)(face)=0;//no slip velocity condition
                        }
                }
                else
                {
                TV_INT cell=face+exterior_cell_offset;
                projection.elliptic_solver->psi_D(cell)=true;
                projection.p(cell)=0;
                }
            }
        }
        else
        {
                for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next())
                        {
                        TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                        projection.elliptic_solver->psi_D(cell)=true;
                        projection.p(cell)=0;
                        }
        }
    }
    for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
    {
     /*  if(time<=3 && source.Lazy_Inside(iterator.Location()))
        {
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==1) face_velocities(iterator.Full_Index())=1;
            else face_velocities(iterator.Full_Index())=0;
        }*/
    }
    for(typename GRID<TV>::CELL_ITERATOR iterator(projection.p_grid);iterator.Valid();iterator.Next())
    {
        if(levelset.phi(iterator.Cell_Index())>0)
        {
        projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;
        projection.p(iterator.Cell_Index())=0;
        }
    }
    if(projection.elliptic_solver->mpi_grid){
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(projection.elliptic_solver->psi_D,1,false);
        projection.elliptic_solver->mpi_grid->Exchange_Boundary_Cell_Data(projection.p,1,false);}
}

//#####################################################################
// 
//#####################################################################



template<class TV>
void LAYOUT<TV>::Write_Output(const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    std::string output_directory="output"; 
    FILE_UTILITIES::Create_Directory("output/"+STRING_UTILITIES::Value_To_String(frame));
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",levelset);
    //FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",mac_grid);
    //FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",levelset);
    if(write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);}

}
template<class TV> void LAYOUT<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/levelset",levelset);
    std::string filename;
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading pressure "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,projection.p);}
}

//
template class LAYOUT<VECTOR<float,2> >;
template class LAYOUT<VECTOR<float,3> >;
