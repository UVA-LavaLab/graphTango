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
	double totTime = 0.0;
	while (!file.eof()) {
		readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);

		t.Start();
		ds->update(el);
		t.Stop();

		totTime += t.Seconds();
		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;

		batch_id++;
	}
	file.close();

	ofstream logFile("time.csv");
	logFile << totTime << endl;
	cout << "Total insertion time: " << totTime << endl;

	t.Start();
	alg.performAlg();
	t.Stop();

	totTime = t.Seconds();
	logFile << totTime << endl;
	cout << "Alg time: " << totTime << endl;

	stringstream ss;
	ss << opts.filename << ".del";
	file.open(ss.str());
	if (!file.is_open()) {
		cout << "Couldn't open file " << ss.str() << endl;
		exit(-1);
	}

	batch_id = 0;
	totTime = 0.0;
	while (!file.eof()) {
		readBatchFromCSV(el, file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);

		t.Start();
		ds->update(el);
		t.Stop();

		totTime += t.Seconds();
		cout << "Deleted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;
		batch_id++;
	}
	file.close();

	logFile << totTime << endl;
	logFile.close();
	cout << "Total deletion time: " << totTime << endl;

	ds->print();
#ifdef CALC_EDGE_TOUCHED
	cout << "EDGES TOUCHED: " << g_edge_touched << endl;
#endif
#ifdef CALC_TYPE_SWITCH
	cout << "Switch count: " << ds->switchCnt << endl;
#endif
}

