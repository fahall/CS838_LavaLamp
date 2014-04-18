#pragma once // compiler directive: only open this file once -- shortens compile times

#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV>
class LAYOUT 
{
public:
	typedef typename TV::SCALAR T;

    const STREAM_TYPE stream_type;
	const int number_of_frames; 	// Total number of frames
    const T frame_time;          	// Frame (snapshot) interval

    LAYOUT(const STREAM_TYPE stream_type_input);
	~LAYOUT();
    void Initialize();
    T Maximum_Dt() const;
    void Write_Output(const int frame) const;
};


}
