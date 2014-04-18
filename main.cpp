#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>

#include "LAYOUT.h"
#include "DRIVER.h"


using namespace PhysBAM;



int main(int argc, char* argv[]) {

    typedef float T;
    typedef float RW;
    typedef VECTOR<T,2> TV;

    RW rw=RW();
    STREAM_TYPE stream_type(rw);



    LOG::Initialize_Logging();

    LAYOUT<TV> layout(stream_type);
    DRIVER<TV> driver(layout);
    driver.Run();

    LOG::Finish_Logging();
}
