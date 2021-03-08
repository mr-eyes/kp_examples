#include "algorithms.hpp"

int main(int argc, char ** argv){


    int kSize=31;
    int chunkSize = 100;
    // string genes_file = "/home/mabuelanin/dib-dev/kp_examples/indexing/data/gencode.v37.pc_transcripts.fa";
    // string genes_file = "/home/mabuelanin/dib-dev/kp_examples/indexing/data/gencode.v37.pc_transcripts.fa.split/gencode.v37.pc_transcripts.part_001.fa";
    string genes_file = argv[1];
    kDataFrame * genesFrame = new kDataFrameMQF(kSize);
    kmerDecoder * KMERS = kProcessor::initialize_kmerDecoder(genes_file, chunkSize, "kmers", {{"k_size", kSize}});

    cerr << "Indexing started ..." << endl;
    // kProcessor::index(KMERS, genes_file+".names", genesFrame);
    kProcessor::index(genesFrame,{{"k_size", kSize}}, genes_file, chunkSize, genes_file+".names");


    genesFrame->save("large_test");
    cerr << "Indexing done ..." << endl;


    return 0;
}