#include "LAYOUT.h"

using namespace PhysBAM;

template<class TV>
LAYOUT<TV>::LAYOUT(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input), number_of_frames(200),frame_time(.05)
{}

template<class TV>
LAYOUT<TV>::~LAYOUT()
{}


template<class TV>
void LAYOUT<TV>::Initialize()
{}


template<class TV>
typename TV::SCALAR LAYOUT<TV>::Maximum_Dt() const
{
    return 1e-3;
}


template<class TV>
void LAYOUT<TV>::Write_Output(const int frame) const
{
    FILE_UTILITIES::Create_Directory("output/"+STRING_UTILITIES::Value_To_String(frame));
    // collection.Write(stream_type,"output",frame,0,true);
}


template class LAYOUT<VECTOR<float,2> >;
template class LAYOUT<VECTOR<float,3> >;
