#include <gtest/gtest.h>
#include "../src/common.h"
#include "../src/filter.h"

TEST(FilerTest, trimAndCut) {
    Read r("@name",
            "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT",
            "+",
            "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E");
    Options opt;
    opt.qualityCut.enabledFront = true;
    opt.qualityCut.enabledTail = true;
    opt.qualityCut.windowSizeFront = 4;
    opt.qualityCut.qualityFront = 20;
    opt.qualityCut.windowSizeTail = 4;
    opt.qualityCut.qualityTail = 20;
    Filter filter(&opt);
    int frontTrimmed = 0;
    Read* ret = filter.trimAndCut(&r, 0, 1, frontTrimmed);

    EXPECT_EQ(*ret->mSeq, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    EXPECT_EQ(*ret->mQuality, "CCCCCCCCCCC////CCCCCCCCCCCCC");
}

TEST(FilerTest, gcContentFilter) {
    Options opt;
    opt.qualfilter.enabled = false;
    opt.lengthFilter.enabled = false;
    opt.gcContentFilter.enabled = true;
    opt.gcContentFilter.min = 25.0;
    opt.gcContentFilter.max = 75.0;

    Filter filter(&opt);

    Read lowGc("@low_gc",
               "AAAAAAAAAATTTTTTTTTT",
               "+",
               "IIIIIIIIIIIIIIIIIIII");
    EXPECT_EQ(filter.passFilter(&lowGc), FAIL_GC_CONTENT);

    Read highGc("@high_gc",
                "GGGGGGGGGGCCCCCCCCCC",
                "+",
                "IIIIIIIIIIIIIIIIIIII");
    EXPECT_EQ(filter.passFilter(&highGc), FAIL_GC_CONTENT);

    Read passing("@passing_gc",
                 "AAAAAGGGGGCCCCCTTTTT",
                 "+",
                 "IIIIIIIIIIIIIIIIIIII");
    EXPECT_EQ(filter.passFilter(&passing), PASS_FILTER);
}

TEST(FilerTest, bestReadSegment) {
    Options opt;
    Filter filter(&opt);

    Read read("@segment",
              "AAAACCCCGGGGTTTT",
              "+",
              "((((IIIIIIII((((");

    int trimmedBases = 0;
    Read* segment = filter.keepBestReadSegment(&read, 4, 20, trimmedBases);
    ASSERT_NE(segment, nullptr);
    EXPECT_NE(segment, &read);
    EXPECT_EQ(*segment->mName, "@best-segment-segment");
    EXPECT_EQ(*segment->mSeq, "CCCCGGGG");
    EXPECT_EQ(*segment->mQuality, "IIIIIIII");
    EXPECT_EQ(trimmedBases, 8);
    delete segment;
}

TEST(FilerTest, bestReadSegmentNoChange) {
    Options opt;
    Filter filter(&opt);

    Read read("@segment",
              "AAAACCCCGGGG",
              "+",
              "IIIIIIIIIIII");

    int trimmedBases = 0;
    Read* segment = filter.keepBestReadSegment(&read, 4, 20, trimmedBases);
    EXPECT_EQ(segment, &read);
    EXPECT_EQ(trimmedBases, 0);
}
