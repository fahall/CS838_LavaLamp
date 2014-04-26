#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>

#include "LAYOUT.h"
#include "DRIVER.h"


using namespace PhysBAM;



int main(int argc, char* argv[]) {

    typedef float T;
    typedef float RW;
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::dimension> TV_INT;
    PARSE_ARGS parse_args;
    parse_args.Add_Integer_Argument("-restart",0,"restart frame");
    parse_args.Add_Integer_Argument("-scale",100,"fine scale grid resolution");
    parse_args.Add_Integer_Argument("-substep",-1,"output-substep level");
    parse_args.Add_Integer_Argument("-e",1,"last frame");
    parse_args.Add_Integer_Argument("-threads",1,"number of threads");
    parse_args.Add_Option_Argument("-3d","run in 3 dimensions");

    RW rw=RW();
    STREAM_TYPE stream_type(rw);
    //TV = PhysBAM::VECTOR<float, 2>;
//PhysBAM::DRIVER<TV>::DRIVER(PhysBAM::LAYOUT<TV>& layout)

    LAYOUT<TV> *layout=new LAYOUT<TV>(stream_type);


    LOG::Initialize_Logging();

    //LAYOUT<TV> *layout(stream_type);
    DRIVER<TV> driver(*layout);

    int scale=parse_args.Get_Integer_Value("-scale"); //scale initialized
    RANGE<TV> range(TV(),TV::All_Ones_Vector()); // range intialized

    TV_INT counts=TV_INT::All_Ones_Vector()*scale;
   // WATER_DRIVER<TV> driver(*layout);
    layout->Initialize_Grid(counts,range);
    driver.Run();
    TV point1=TV::All_Ones_Vector()*(.65),point2=TV::All_Ones_Vector()*.75;//gives 2 points created
point1(1)=0;//moves x plane
point2(1)=.05;//moves xplane
point1(2)=0;//extends y plane

    layout->source.min_corner=point1;layout->source.max_corner=point2;//the points are now set as upper and lower corner of injector

    LOG::Finish_Logging();

return 0;
}
