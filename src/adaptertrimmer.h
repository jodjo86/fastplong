#ifndef ADAPTER_TRIMMER_H
#define ADAPTER_TRIMMER_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "filterresult.h"
#include "options.h"
#include "read.h"

using namespace std;

class AdapterTrimmer{
public:
    AdapterTrimmer();
    ~AdapterTrimmer();

    static int trimBySequenceStart(Read* r1, FilterResult* fr, string& adapter, double edMax = 0.3, int trimmingExtension = 10);
    static int trimBySequenceEnd(Read* r1, FilterResult* fr, string& adapter, double edMax = 0.3, int trimmingExtension = 10);
    static int trimByMultiSequences(Read* r1, FilterResult* fr, vector<string>& adapterList, double edMax = 0.3, int trimmingExtension = 10);
    static bool findMiddleAdapters(Read* r, string& startAdater, string& endAdapter, int& start, int& len, double edMax = 0.3, int trimmingExtension = 10);
    static bool findNextMiddleAdapter(Read* r, vector<string>& adapters, int searchStart, int& adapterStart, int& adapterLength, double edMax = 0.3);
    static vector<Read*> splitByMiddleAdapters(Read* r, FilterResult* fr, vector<string>& adapters, double edMax = 0.3, int trimmingExtension = 10, int minSegmentLength = 20, int* removedBases = NULL);
    static int searchAdapter(string* read, const string& adapter, double edMax = 0.3, int searchStart = 0, int searchLen = -1, bool asLeftAsPossible = false, bool asRightAsPossible = false);
};


#endif
