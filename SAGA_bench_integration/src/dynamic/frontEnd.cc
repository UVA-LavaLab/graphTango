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

	int batch_id = 0;
	NodeID lastAssignedNodeID = -1;
	MapTable VMAP;
	ofstream updF("Update.csv");

	Timer t;
	dataStruc *ds = createDataStruc(opts.type, opts.weighted, opts.directed, opts.num_nodes, opts.num_threads);
	Algorithm alg(opts.algorithm, ds, opts.type);

//	while (!file.eof()) {
//		EdgeList el = readBatchFromCSV(file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);
//
//		t.Start();
//		ds->update(el);
//		t.Stop();
//
//		updF << t.Seconds() << endl;
//		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;
//
//		alg.performAlg();
//
//		batch_id++;
//	}
//	updF.close();


	while (!file.eof()) {
		EdgeList el = readBatchFromCSV(file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);
		ds->update(el);
		cout << "Inserted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;
		batch_id++;
	}

	file.close();

////	stringstream ss;
////	ss << opts.filename << ".del";
////	file.open(ss.str());
////	if (!file.is_open()) {
////		cout << "Couldn't open file " << ss.str() << endl;
////		exit(-1);
////	}
////
////	batch_id = 0;
////	while (!file.eof()) {
////		EdgeList el = readBatchFromCSV(file, opts.batch_size, batch_id, opts.weighted, VMAP, lastAssignedNodeID);
////
////		t.Start();
////		ds->update(el);
////		t.Stop();
////
////		updF << t.Seconds() << endl;
////		cout << "Deleted Batch " << batch_id << ": Nodes " << ds->num_nodes << ", Edges " << ds->num_edges << endl;
////
////		alg.performAlg();
////
////		batch_id++;
////	}
////	updF.close();



	ds->print();
}
