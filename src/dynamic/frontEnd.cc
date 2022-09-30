#include <unistd.h>
#include <fstream>
#include <cstring>
#include <mutex>
#include <thread>
#include <pthread.h>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "topAlg.h"
#include "topDataStruc.h"
#include "fileReader.h"
#include "parser.h"
#include "../common/timer.h"


using namespace std;
/* Main thread that launches everything else */

ofstream algF("Alg.csv");

typedef struct {
	i32 src;
	i32 dst;
} BinEdge;

u64 readBatchFromBin(u64 numNodes, u64 numEdges, BinEdge* buffer, EdgeList& el, FILE* in, int batchSize, bool weighted){
	el.clear();
	static u64 readEdges = 0;
	static i32 readNodeId = -1;

	u64 remEdges = numEdges - readEdges;
	if(remEdges > batchSize){
		remEdges = batchSize;
	}

	assert(fread(buffer, sizeof(BinEdge), remEdges, in) == remEdges);

	for(u64 i = 0; i < remEdges; i++){
		Edge e;
		e.source = buffer[i].src;
		e.destination = buffer[i].dst;
		if(e.source < 0){
			//deletion
			assert(e.destination < 0);
			e.source = e.source ^ 0xffffffff;
			e.destination = e.destination ^ 0xffffffff;
			e.isDelete = true;
		}
		if(weighted){
			e.weight = rand() % 256;
		}

		//update max nodes so far
		if(e.source > readNodeId){
			readNodeId = e.source;
			e.sourceExists = false;
		}
		if(e.destination > readNodeId){
			readNodeId = e.destination;
			e.destExists = false;
		}
		el.push_back(e);
	}

	assert(readNodeId < numNodes);
	readEdges += remEdges;
	return readNodeId + 1;
}

int main(int argc, char *argv[]) {
	ios_base::sync_with_stdio(false);

	cmd_args opts = parse(argc, argv);

	FILE* in = fopen(opts.filename.c_str(), "r");
	if(!in) {
		perror("Cannot open file.");
		exit(-1);
	}

	u64 numNodes;
	u64 numEdges;
	u64 sz;
	sz = fread(&numNodes, 1, 8, in);
	assert(sz == 8);
	sz = fread(&numEdges, 1, 8, in);
	assert(sz == 8);
	srand(42);

	cout << "[" << opts.filename << "] Nodes: " << numNodes << endl;
	cout << "[" << opts.filename << "] Edges: " << numEdges << endl;

	int batch_id = 0;
	EdgeList el;
	el.reserve(opts.batch_size);

	Timer t;
	dataStruc *ds = createDataStruc(opts.type, opts.weighted, opts.directed, numNodes, opts.num_threads);
	Algorithm alg(opts.algorithm, ds, opts.type);

	ofstream updF("Update.csv");

	BinEdge* buffer = (BinEdge*)malloc(opts.batch_size * sizeof(BinEdge));
	const u64 numBatch = ceil(numEdges * 1.0 / opts.batch_size);
	for(u64 batch_id = 0; batch_id < numBatch; batch_id++){
		readBatchFromBin(numNodes, numEdges, buffer, el, in, opts.batch_size, opts.weighted);
		t.Start();
		ds->update(el);
		t.Stop();

		updF << t.Seconds() << endl;
		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;

		//alg.performAlg();
	}
	free(buffer);
	updF.close();


//	while (!file.eof()) {
//		readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);
//		ds->update(el);
//		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;
//		//cout << "ins," << ((ds->num_edges * 1.0) / ds->num_nodes) << endl;
//		batch_id++;
//	}
//
//	file.close();

//	stringstream ss;
//	ss << opts.filename << ".del";
//	file.open(ss.str());
//	if (!file.is_open()) {
//		cout << "Couldn't open file " << ss.str() << endl;
//		exit(-1);
//	}
//
//	batch_id = 0;
//	while (!file.eof()) {
//		readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);
//
//		t.Start();
//		ds->update(el);
//		t.Stop();
//
//		updF << t.Seconds() << endl;
//		//cout << "del," << ((ds->num_edges * 1.0) / ds->num_nodes) << endl;
//		cout << "Deleted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;
//
//		alg.performAlg();
//
//		batch_id++;
//	}
//	updF.close();

	ds->print();
#ifdef CALC_EDGE_TOUCHED
	cout << "EDGES TOUCHED: " << g_edge_touched << endl;
#endif
#ifdef CALC_TYPE_SWITCH
	cout << "Switch count: " << ds->switchCnt << endl;
#endif
}
