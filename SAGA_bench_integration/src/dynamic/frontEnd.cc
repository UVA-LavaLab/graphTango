#include <unistd.h>
#include <fstream>
#include <cstring>
#include <mutex>
#include <thread>
#include <pthread.h>
#include <sstream>

#include "topAlg.h"
#include "topDataStruc.h"
//#include "builder.h"
#include "fileReader.h"
#include "parser.h"
#include "../common/timer.h"


using namespace std;
/* Main thread that launches everything else */

ofstream algF("Alg.csv");

int main(int argc, char *argv[]) {
	cmd_args opts = parse(argc, argv);
	ifstream file(opts.filename);
	if (!file.is_open()) {
		cout << "Couldn't open file " << opts.filename << endl;
		exit(-1);
	}

	LIKWID_MARKER_INIT;

	#pragma omp parallel
	{
	  LIKWID_MARKER_THREADINIT;
	}

	int batch_id = 0;
	NodeID lastAssignedNodeID = -1;
	MapTable VMAP;
	EdgeList el;
	el.reserve(opts.batch_size);

	Timer t;
	dataStruc *ds = createDataStruc(opts.type, opts.weighted, opts.directed, opts.num_nodes, opts.num_threads);
	Algorithm alg(opts.algorithm, ds, opts.type);


	ofstream updF("Update.csv");

	u32 maxIn = 0;
	u32 maxOut = 0;

	unordered_map<u32, u32> inMap;
	unordered_map<u32, u32> outMap;

	while (!file.eof()) {
		readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);

		t.Start();
		ds->update(el);
		t.Stop();

		for(const auto it : el){
			u32 src = it.source;
			u32 dst = it.destination;
			inMap[dst]++;
			outMap[src]++;
			if(inMap[dst] > maxIn){
				maxIn = inMap[dst];
			}
			if(outMap[src] > maxOut){
				maxOut = outMap[src];
			}
		}
		inMap.clear();
		outMap.clear();

		updF << t.Seconds() << endl;
		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;

		alg.performAlg();

		batch_id++;
	}
	updF.close();


	cout << "Max in: " << maxIn << endl;
	cout << "Max Out: " << maxOut << endl;

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
	LIKWID_MARKER_CLOSE;
}
