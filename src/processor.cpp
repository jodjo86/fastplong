#include "processor.h"
#include "seprocessor.h"

Processor::Processor(Options* opt){
    mOptions = opt;
}


Processor::~Processor(){
}

bool Processor::process(ProcessingResult* result) {
    SingleEndProcessor p(mOptions);
    p.process(result);

    return true;
}
