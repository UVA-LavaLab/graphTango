#ifndef FILEREADER_H_
#define FILEREADER_H_

#include <sstream>
#include <fstream>
#include <cstdio>

#include "types.h"

using namespace std;

/* 
1) Read edge streams from files
2) Maintain MapTable
3) Assign logical IDs
4) Assign batch IDs 

TO DO: Add support for other file types
*/

Edge convertCSVLineIntoEdge(const char delim, const string& line, bool weighted){
    // extract the numbers from the line-> a,b,c,d
    stringstream ss(line); // put the line in an internal stream     
    Edge e;
    string data;
    
    getline(ss, data, delim); e.source = stol(data);
    getline(ss, data, delim); e.destination = stol(data);
    //getline(ss, data, delim); /*time = stol(data);*/

    if(weighted){
        //getline(ss, data, delim);
        //e.weight = stol(data);
	e.weight = (rand() % 8) + 8;
    }

    if(line[0] == '-'){
    	e.source = e.source * (-1);
    	e.isDelete = true;
    }
    return e;    
}

// return true if a mapping is found, otherwise false
bool assignLogicalID(NodeID& n, MapTable& VMap, NodeID& lastAssignedLogicalID){
    if(VMap.empty()){
        // this is the first vertex ever         
        VMap[n] = 0;  
        n = 0;
        lastAssignedLogicalID=0;
        return false;        
    }
    
    MapTable::iterator it;
    it = VMap.find(n);
    if(it != VMap.end()){
        // vertex exists
        n = it->second;
        return true;
    }
    
    else{
        // vertex does not exist   
        lastAssignedLogicalID++;           
        VMap[n] = lastAssignedLogicalID; 
        n = lastAssignedLogicalID;          
        return false;
    }
}


void readBatchFromCSV(EdgeList& el, ifstream& in, int batchSize, int batch_id, bool weighted, MapTable& VMap, NodeID& lastAssignedLogicalID){
	el.clear();
    int edgecount = 0;
    string line;      

    while(getline(in, line)){
        if(line != ""){          
            Edge e = convertCSVLineIntoEdge(' ', line, weighted);
            if(assignLogicalID(e.source, VMap, lastAssignedLogicalID)) e.sourceExists = true;
            if(assignLogicalID(e.destination, VMap, lastAssignedLogicalID)) e.destExists = true;
            //e.batch_id = batch_id;
            //e.lastAssignedId = lastAssignedLogicalID;
            el.push_back(e);
            edgecount++;  
            if(edgecount == batchSize) break;   
        }             
    }
}

#endif  // FILEREADER_H_
