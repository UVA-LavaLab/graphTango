#include <unistd.h>
#include <fstream>
#include <cstring>
#include <mutex>
#include <thread>
#include <pthread.h>
#include <sstream>

#include "topAlg.h"
#include "topDataStruc.h"
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

	int batch_id = 0;
	NodeID lastAssignedNodeID = -1;
	MapTable VMAP;
	EdgeList el;
	el.reserve(opts.batch_size);

	Timer t;
	dataStruc *ds = createDataStruc(opts.type, opts.weighted, opts.directed, opts.num_nodes, opts.num_threads);
	Algorithm alg(opts.algorithm, ds, opts.type); 
	std::string outputFile;
        if (!opts.outFileName.empty())
		outputFile = opts.outFileName;
	else
		outputFile = "update.csv";
	ofstream updF(opts.outFileName);
	while (!file.eof()) {
		readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);

		t.Start();
		ds->update(el);
		t.Stop();

		updF << t.Seconds() << endl;
		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;


		batch_id++;
	}
	updF.close();

#ifndef ENABLE_PROFILING /*No need to run algorithm if profiling non algorithm specific things, e.g., type mapping */
		alg.performAlg();
#endif


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
	if(ds){
		delete ds;
	}

#ifdef CALC_EDGE_TOUCHED
	cout << "EDGES TOUCHED: " << g_edge_touched << endl;
#endif

#ifdef CALC_PROBING_DISTANCE
	cout << "Average probing distance: " << 1.0 * g_num_probes / g_num_type3 << endl;
#endif
}

