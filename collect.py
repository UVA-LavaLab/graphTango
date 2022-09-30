#!/usr/bin/env python3

import sys
import os

algoMap = {
	"bfsfromscratch" 			: "bfs_fs",
	"bfsdyn"					: "bfs",
	"prfromscratch"				: "pr_fs",
	"prdyn" 					: "pr",
	"ssspfromscratch" 			: "sssp_fs",
	"ssspdyn" 					: "sssp",
	"sswpfromscratch"			: "sswp_fs",
	"sswpdyn"					: "sswp",
	"mcfromscratch"				: "mc_fs",
	"mcdyn"						: "mc",
	"ccfromscratch"				: "cc_fs",
	"ccdyn"						: "cc"
}

datasetMap = {
	"LiveJournal.csv"			: "LiveJournal",
	"WikiTalk.csv"				: "Talk",
	"orkut.csv"					: "Orkut",
	"wiki-topcats.csv"			: "Wiki",
	"rmat_1_1.csv"				: "1",
	"rmat_1_2.csv"				: "2",
	"rmat_1_4.csv"				: "4",
	"rmat_1_8.csv"				: "8",
	"rmat_1_16.csv"				: "16",
	"rmat_1_32.csv"				: "32",
	"rmat_1_64.csv"				: "64",
	"rmat_1_128.csv"				: "128",
	"rmat_1_256.csv"				: "256",
}

dsMap = {
	#"stinger"					: "Stinger",
	#"adListShared"				: "AdListShared",
	#"adListChunked"				: "AdListChunked",
	#"degAwareRHH"				: "DegAwareRHH",
	#"graphite"					: "Ours",
	
	#"graphite_2"				: "Threshold 2",
	#"graphite_4"				: "Threshold 4",
	#"graphite_8"				: "Threshold 8",
	#"graphite_16"				: "Threshold 16",
	#"graphite_32"				: "Threshold 32",
	#"graphite_64"				: "Threshold 64",
	#"graphite_128"				: "Threshold 128",
	#"graphite_256"				: "Hybrid",
	#"graphite_512"				: "Threshold 512",
	#"graphite_1024"				: "Threshold 1024",
	#"graphite_hash"				: "Hash only",
	#"graphite_linear"			: "Linear only",
	
	#"graphite_512_1"			: "Th_512_Min_1",
	#"graphite_512_2"			: "Th_512_Min_2",
	#"graphite_512_4"			: "Th_512_Min_4",
	#"graphite_512_8"			: "Th_512_Min_8",
	#"graphite_512_16"			: "Th_512_Min_16",
	#"graphite_512_32"			: "Th_512_Min_32",
	
	#"graphite_128_grp"			: "Threshold 128 grp",
	#"graphite_256_grp"			: "Threshold 256 grp",
	#"graphite_512_grp"			: "Threshold 512 grp",
	#"graphite_1024_grp"			: "Threshold 1024 grp",
	
	#"graphite_512_lock"			: "Threshold 512 lock",
	#"graphite_512_tight" 		: "Threshold 512 tight",
	#"graphite_512_cfh" 			: "Threshold 512 cfh",
	#"graphite_32_cfh" 			: "Threshold 32 cfh",
	#"graphite_cfh_only" 		: "cfh only",
	#"graphite_cfh_2cache"		: "cfh 2 lines",
	#"graphite_chf_2"			: "CFH",

	#"graphite_cfh_256"			: "CFH + adjList (256)",
	#"graphite_cfh_128"			: "CFH + adjList (128)",
	#"graphite_cfh_64"			: "CFH + adjList",
	#"graphite_cfh_32"			: "CFH + adjList (32)",
	#"graphite_cfh_16"			: "CFH + adjList (16)",
	
	#"del_graphite_cfh_64"		: "CFH + adjList",
	#"del_stinger"				: "Stinger",
	#"del_adListShared"			: "AdListShared",
	#"del_adListChunked"			: "AdListChunked",
	#"del_degAwareRHH"			: "DegAwareRHH",

	#"threads2_graphite"			: "GraphTango 2Th",
	#"threads2_stinger"			: "Stinger 2Th",
	#"threads2_adListShared"		: "AdListShared 2Th",
	#"threads2_adListChunked"	: "AdListChunked 2Th",
	
	#"threads4_graphite"			: "GraphTango 4Th",
	#"threads4_stinger"			: "Stinger 4Th",
	#"threads4_adListShared"		: "AdListShared 4Th",
	#"threads4_adListChunked"	: "AdListChunked 4Th",
	#"threads4_degAwareRHH"		: "DegAwareRHH 4Th",	
	
	#"threads8_graphite"			: "GraphTango 8Th",
	#"threads8_stinger"			: "Stinger 8Th",
	#"threads8_adListShared"		: "AdListShared 8Th",
	#"threads8_adListChunked"	: "AdListChunked 8Th",
	#"threads8_degAwareRHH"		: "DegAwareRHH 8Th",
	
	
	#-------------------TH1 sweep for GTBalanced--------------
	#"gtbalanced_th1_8"		: "8",
	#"gtbalanced_th1_16"		: "16",
	#"gtbalanced_th1_32"		: "32",
	#"gtbalanced_th1_64"		: "64",
	#"gtbalanced_th1_128"	: "128",
	#"gtbalanced_th1_256"	: "256",
	#"gtbalanced_th1_512"	: "512",
	

	#-------------------Core sweep for AdjListShared--------------
	#"adListShared_cores_2"			: "AdListShared 2",
	#"adListShared_cores_4"			: "AdListShared 4",
	#"adListShared_cores_6"			: "AdListShared 6",
	#"adListShared_cores_8"			: "AdListShared 8",
	#"adListShared_cores_10"			: "AdListShared 10",
	#"adListShared_cores_12"			: "AdListShared 12",
	

	#-------------------Core sweep for AdjListChunked--------------
	#"adListChunked_cores_2"		: "AdListChunked 2",
	#"adListChunked_cores_4"		: "AdListChunked 4",
	#"adListChunked_cores_6"		: "AdListChunked 6",
	#"adListChunked_cores_8"		: "AdListChunked 8",
	#"adListChunked_cores_10"		: "AdListChunked 10",
	#"adListChunked_cores_12"		: "AdListChunked 12",

	

	#-------------------Core sweep for Stinger--------------
	#"stinger_cores_2"				: "Stinger 2",
	#"stinger_cores_4"				: "Stinger 4",
	#"stinger_cores_6"				: "Stinger 6",
	#"stinger_cores_8"				: "Stinger 8",
	#"stinger_cores_10"				: "Stinger 10",
	#"stinger_cores_12"				: "Stinger 12",

	
	#-------------------Core sweep for DAH--------------
	#"degAwareRHH_cores_2"			: "DAH 2",
	#"degAwareRHH_cores_4"			: "DAH 4",
	#"degAwareRHH_cores_6"			: "DAH 6",
	#"degAwareRHH_cores_8"			: "DAH 8",
	#"degAwareRHH_cores_10"			: "DAH 10",
	#"degAwareRHH_cores_12"			: "DAH 12",

	

	#-------------------Core sweep for GTBalanced--------------
	#"gtbalanced_th1_32_cores_2"		: "TangoBalanced 2",
	#"gtbalanced_th1_32_cores_4"		: "TangoBalanced 4",
	#"gtbalanced_th1_32_cores_6"		: "TangoBalanced 6",
	#"gtbalanced_th1_32_cores_8"		: "TangoBalanced 8",
	#"gtbalanced_th1_32_cores_10"	: "TangoBalanced 10",
	#"gtbalanced_th1_32"	: "TangoBalanced 12",

	
	#-------------------Main comparison---------------
	#"adListShared_cores_12"			: "AdListShared",
	#"adListChunked_cores_12"		: "AdListChunked",
	#"stinger_cores_12"				: "Stinger",
	#"degAwareRHH_cores_12"			: "DAH",
	"gtbalanced_th1_32"				: "GraphTango",
	"graphTango"				: "temp",
		
		
	#-------------------Delete comparison---------------
	#"adListShared_cores_12_del"			: "AdListShared",
	#"adListChunked_cores_12_del"		: "AdListChunked",
	#"stinger_cores_12_del"				: "Stinger",
	#"degAwareRHH_cores_12_del"			: "DAH",
	#"gtbalanced_th1_32_del"	: "GraphTango",
		
	
	#-------------------Comparison of sub-optimal versions---------------
	#"gtbalanced_th1_32_type3"	: "Type3 only",
	#"gtbalanced_th1_32_abseil"	: "Abseil",
	#"gtbalanced_th1_32_RHH"	: "GT RHH",
	#"gtbalanced_th1_32_tsl_rhh"	: "GT TSL",
	"gtbalanced_th1_32_no_switch"	: "GT No Switch",
	"gtbalanced_th1_32_only_growth"	: "GT Only Growth",
	#"gtbalanced_th1_32_malloc_stdmap"	: "Hybrid only",
	#"gtbalanced_th1_32_stdmap"	: "Hybrid + unordered_map",
	#"gtbalanced_th1_32_malloc"	: "Hybrid + malloc",
	#"gtbalanced_th1_32"			: "GraphTango",
	#"degAwareCFH"				: "degAwareCFH",
	#"degAwareCFH_2"				: "degAwareCFH2",
	#"degAwareRHH_cores_12"			: "DAH",
	#"gtbalanced_th1_32_hugepage" : "GraphTango Hugepage",
	#"graphite"					: "GraphTango NN",
	
	
	#--------------Impact of average degree-------------------------
	#"adListShared_cores_12"			: "AdListShared",
	#"adListChunked_cores_12"		: "AdListChunked",
	#"stinger_cores_12"				: "Stinger",
	#"degAwareRHH_cores_12"			: "DAH",
	#"gtbalanced_th1_32"				: "GraphTango",
}

#activeAlgoList = ["bfsdyn", "prdyn", "ssspdyn", "sswpdyn", "ccdyn", "mcdyn"]
#activeAlgoList = ["bfsdyn", "ssspdyn", "sswpdyn", "ccdyn", "mcdyn"]
activeAlgoList = ["bfsdyn"]
activeDatasetList = ["orkut.csv", "LiveJournal.csv", "wiki-topcats.csv", "WikiTalk.csv"]
#activeDatasetList = ["rmat_1_1.csv", "rmat_1_2.csv", "rmat_1_4.csv", "rmat_1_8.csv", "rmat_1_16.csv", "rmat_1_32.csv", "rmat_1_64.csv", "rmat_1_128.csv", "rmat_1_256.csv"]

def build_insert_csv(rootDir):
	with open(rootDir + "/update.csv", "w") as outf:
		for dataset in activeDatasetList:
			headerPrinted = False
			header = datasetMap[dataset] + ",batch,"
			for ds in dsMap:
				fname = rootDir + '/bfsdyn/' + ds + '/' + dataset + '/Update.csv'
				print("Processing " + fname)
				with open(fname, "r") as csvf:
					lineCount = 0
					currOut = datasetMap[dataset] + ',' + dsMap[ds] + ','
					for line in csvf:
						if headerPrinted == False:
							lineCount = lineCount + 1
							header = header + str(lineCount) + ','
						currOut = currOut + line[:-1] + ','
					currOut = currOut + '\n'
					if headerPrinted == False:
						header = header + '\n'
						outf.write(header)
						headerPrinted = True
					outf.write(currOut)
			outf.write(",\n,\n")
			
			

def build_algo_csv(rootDir):
	for algo in activeAlgoList:
		with open(rootDir + "/" + algoMap[algo] + ".csv", "w") as outf:
			for dataset in activeDatasetList:
				headerPrinted = False
				header = datasetMap[dataset] + ",batch,"
				for ds in dsMap:
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Alg.csv'
					print("Processing " + fname)
					with open(fname, "r") as csvf:
						lineCount = 0
						currOut = datasetMap[dataset] + ',' + dsMap[ds] + ','
						for line in csvf:
							if headerPrinted == False:
								lineCount = lineCount + 1
								header = header + str(lineCount) + ','
							currOut = currOut + line[:-1] + ','
						currOut = currOut + '\n'
						if headerPrinted == False:
							header = header + '\n'
							outf.write(header)
							headerPrinted = True
						outf.write(currOut)
				outf.write(",\n,\n")



def build_summed_algo_csv(rootDir):
	with open(rootDir + "/sum_algo.csv", "w") as outf:
		for algo in activeAlgoList:
			outf.write(",");
			for dataset in activeDatasetList:
				outf.write(',' + datasetMap[dataset])
			outf.write('\n')
			
			for ds in dsMap:
				outf.write(algoMap[algo] + ',' + dsMap[ds])
			
				for dataset in activeDatasetList:
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Alg.csv'
					print("Processing " + fname)
					with open(fname, "r") as csvf:
						tot = 0.0
						for line in csvf:
							tot = tot + float(line)
					outf.write(',' + str(tot))
				outf.write('\n')
			outf.write('\n\n\n')


def build_summed_update_csv(rootDir):
	with open(rootDir + "/sum_update.csv", "w") as outf:
		for algo in activeAlgoList:
			outf.write(",");
			for dataset in activeDatasetList:
				outf.write(',' + datasetMap[dataset])
			outf.write('\n')
			
			for ds in dsMap:
				outf.write(algoMap[algo] + ',' + dsMap[ds])
			
				for dataset in activeDatasetList:
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Update.csv'
					print("Processing " + fname)
					with open(fname, "r") as csvf:
						tot = 0.0
						for line in csvf:
							tot = tot + float(line)
					outf.write(',' + str(tot))
				outf.write('\n')
			outf.write('\n\n\n')


def build_summed_tot_csv(rootDir):
	with open(rootDir + "/total.csv", "w") as outf:
		for algo in activeAlgoList:
			outf.write(",");
			for dataset in activeDatasetList:
				outf.write(',' + datasetMap[dataset])
			outf.write('\n')
			
			for ds in dsMap:
				outf.write(algoMap[algo] + ',' + dsMap[ds])
			
				for dataset in activeDatasetList:
					tot = 0.0
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Update.csv'
					print("Processing " + fname)
					with open(fname, "r") as csvf:				
						for line in csvf:
							tot = tot + float(line)
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Alg.csv'
					print("Processing " + fname)
					with open(fname, "r") as csvf:				
						for line in csvf:
							tot = tot + float(line)							
					outf.write(',' + str(tot))
				outf.write('\n')
			outf.write('\n\n\n')


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("usage: ./collect.py <result_root_dir>")
		exit(-1)
	
	rootDir = sys.argv[1]
	
#	build_insert_csv(rootDir)
#	build_algo_csv(rootDir)
	build_summed_algo_csv(rootDir)
	build_summed_update_csv(rootDir)
	build_summed_tot_csv(rootDir)	

	
	
		
		
