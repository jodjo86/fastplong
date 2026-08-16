#include <gtest/gtest.h>
#include "../src/adaptertrimmer.h"

TEST(AdapterTrimmer, trimBySequenceStart) {
    Read r("@name",
        "AGGTGCTGCGCATACTTTTCCACGGGGATACTACTGGGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAAAGGTGCT",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEE///EEEEEEEEEEEEEEEEEEEEE");
    string adapter = "GCGCATACTTTTCCACGGGGATACTACTG";
    int trimmed = AdapterTrimmer::trimBySequenceStart(&r, NULL, adapter, 0.3, 0);
    EXPECT_EQ(*r.mSeq, "GGTGTTACCGTGGGAATGAATCCTTTTAACCTTAGCAATACGTAAAGGTGCT");

    Read r2("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAAGCGCATACTTTTCCACGGGGA",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEET");
    trimmed = AdapterTrimmer::trimBySequenceEnd(&r2, NULL, adapter, 0.3, 0);
    EXPECT_EQ(*r2.mSeq, "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA");

    // Read read("@name",
    //     "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC",
    //     "+",
    //     "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    // vector<string> adapterList;
    // adapterList.push_back("GCTAGCTAGCTAGCTA");
    // adapterList.push_back("AAATTTCCCGGGAAATTTCCCGGG");
    // adapterList.push_back("ATCGATCGATCGATCG");
    // adapterList.push_back("AATTCCGGAATTCCGG");
    // trimmed = AdapterTrimmer::trimByMultiSequences(&read, NULL, adapterList);
    // if (*read.mSeq != "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG") {
    //     cerr << read.mSeq << endl;
    //     return false;
    // }
}


TEST(AdapterTrimmer, searchAdapterLeft) {
    Read read("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    string adapter = "TTTT";
    int pos = AdapterTrimmer::searchAdapter(
        read.mSeq, adapter, 0.3, 0, -1, true, false);
    EXPECT_EQ(pos, 0);
}

TEST(AdapterTrimmer, searchAdapterLeft2) {
    Read read("@name",
        "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC",
        "+",
        "///EEEEEEEEEEEEEEEEEEEEEEEEEE////EEEEEEEEEEEEE////E////EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    string adapter = "AACC";
    int pos = AdapterTrimmer::searchAdapter(
        read.mSeq, adapter, 0.3, 0, -1, true, false);
    EXPECT_EQ(pos, 4);
}

TEST(AdapterTrimmer, splitByMiddleAdapters) {
    Read read("@name",
        "AAAACCCCGGGGTTTTCCCCGGGGAAAA",
        "+",
        "EEEEEEEEEEEEEEEEEEEEEEEEEEEE");
    vector<string> adapters;
    adapters.push_back("CCCCGGGG");
    int removedBases = 0;

    vector<Read*> segments = AdapterTrimmer::splitByMiddleAdapters(&read, NULL, adapters, 0.0, 0, 4, &removedBases);

    ASSERT_EQ(segments.size(), 3u);
    EXPECT_EQ(*segments[0]->mSeq, "AAAA");
    EXPECT_EQ(*segments[1]->mSeq, "TTTT");
    EXPECT_EQ(*segments[2]->mSeq, "AAAA");
    EXPECT_EQ(*segments[0]->mName, "@chimera-segment-1-name");
    EXPECT_EQ(*segments[1]->mName, "@chimera-segment-2-name");
    EXPECT_EQ(*segments[2]->mName, "@chimera-segment-3-name");
    EXPECT_EQ(removedBases, 16);

    for(size_t i=0; i<segments.size(); i++)
        delete segments[i];
}

TEST(AdapterTrimmer, splitByMiddleAdaptersNoMatchReturnsOriginal) {
    Read read("@name",
        "AAAATTTTAAAA",
        "+",
        "EEEEEEEEEEEE");
    vector<string> adapters;
    adapters.push_back("CCCCGGGG");

    vector<Read*> segments = AdapterTrimmer::splitByMiddleAdapters(&read, NULL, adapters, 0.0, 0, 4);

    ASSERT_EQ(segments.size(), 1u);
    EXPECT_EQ(segments[0], &read);
}
