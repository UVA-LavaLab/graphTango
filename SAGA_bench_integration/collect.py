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
	"wiki-topcats.csv"			: "Wiki"
}

dsMap = {
	"stinger"					: "Stinger",
	"adListShared"				: "AdListShared",
	"adListChunked"				: "AdListChunked",
	"degAwareRHH"				: "DegAwareRHH",
	"graphite"					: "Ours",
	
	"graphite_2"				: "Threshold 2",
	"graphite_4"				: "Threshold 4",
	"graphite_8"				: "Threshold 8",
	"graphite_16"				: "Threshold 16",
	"graphite_32"				: "Threshold 32",
	"graphite_64"				: "Threshold 64",
	"graphite_128"				: "Threshold 128",
	"graphite_256"				: "Hybrid",
	"graphite_512"				: "Threshold 512",
	"graphite_1024"				: "Threshold 1024",
	"graphite_hash"				: "Hash only",
	"graphite_linear"			: "Linear only",
	
	"graphite_512_1"			: "Th_512_Min_1",
	"graphite_512_2"			: "Th_512_Min_2",
	"graphite_512_4"			: "Th_512_Min_4",
	"graphite_512_8"			: "Th_512_Min_8",
	"graphite_512_16"			: "Th_512_Min_16",
	"graphite_512_32"			: "Th_512_Min_32",
	
	"graphite_128_grp"			: "Threshold 128 grp",
	"graphite_256_grp"			: "Threshold 256 grp",
	"graphite_512_grp"			: "Threshold 512 grp",
	"graphite_1024_grp"			: "Threshold 1024 grp",
	
	"graphite_512_lock"			: "Threshold 512 lock",
	"graphite_512_tight" 		: "Threshold 512 tight",
	"graphite_512_cfh" 			: "Threshold 512 cfh",
	"graphite_32_cfh" 			: "Threshold 32 cfh",
	"graphite_cfh_only" 		: "cfh only",
	"graphite_cfh_2cache"		: "cfh 2 lines",
	"graphite_chf_2"			: "CFH",

	"graphite_cfh_256"			: "CFH + adjList (256)",
	"graphite_cfh_128"			: "CFH + adjList (128)",
	"graphite_cfh_64"			: "CFH + adjList",
	"graphite_cfh_32"			: "CFH + adjList (32)",
	"graphite_cfh_16"			: "CFH + adjList (16)",
	
	"del_graphite_cfh_64"		: "CFH + adjList",
	"del_stinger"				: "Stinger",
	"del_adListShared"			: "AdListShared",
	"del_adListChunked"			: "AdListChunked",
	"del_degAwareRHH"			: "DegAwareRHH",

	"threads2_graphite"			: "GraphTango 2Th",
	"threads2_stinger"			: "Stinger 2Th",
	"threads2_adListShared"		: "AdListShared 2Th",
	"threads2_adListChunked"	: "AdListChunked 2Th",
	
	"threads4_graphite"			: "GraphTango 4Th",
	"threads4_stinger"			: "Stinger 4Th",
	"threads4_adListShared"		: "AdListShared 4Th",
	"threads4_adListChunked"	: "AdListChunked 4Th",
	"threads4_degAwareRHH"		: "DegAwareRHH 4Th",	
	
	"threads8_graphite"			: "GraphTango 8Th",
	"threads8_stinger"			: "Stinger 8Th",
	"threads8_adListShared"		: "AdListShared 8Th",
	"threads8_adListChunked"	: "AdListChunked 8Th",
	"threads8_degAwareRHH"		: "DegAwareRHH 8Th",		
}

#activeAlgoList = ["bfsdyn", "prdyn", "ssspdyn", "sswpdyn", "ccdyn", "mcdyn"]
activeAlgoList = ["bfsdyn"]
activeDatasetList = ["orkut.csv", "LiveJournal.csv", "wiki-topcats.csv", "WikiTalk.csv"]
#activeDSList = ["adListShared", "adListChunked", "stinger", "degAwareRHH", "graphite"]
activeDSList = [
	"degAwareRHH",
	"stinger",
	"adListShared",
	"adListChunked",
	#"graphite_2",
	#"graphite_4",
	#"graphite_8",
	#"graphite_16",
	#"graphite_32",
	#"graphite_64",
	#"graphite_128",
	#"graphite_256",
	#"graphite_512",
	#"graphite_1024",
	#"graphite_hash",
	#"graphite_linear",
	#"graphite_512_1",
	#"graphite_512_2",
	#"graphite_512_4",
	#"graphite_512_8",
	#"graphite_512_16",
	#"graphite_512_32",
	#"graphite_128_grp",
	#"graphite_256_grp",
	#"graphite_512_grp",
	#"graphite_1024_grp",
	#"graphite_512_lock",
	#"graphite_512_tight",
	#"graphite_512_cfh",
	#"graphite_32_cfh",
	#"graphite_cfh_only",
	#"graphite_cfh_2cache",
	#"graphite_chf_2",
	"graphite_cfh_64",
	
	#"del_graphite_cfh_64",
	#"del_stinger",
	#"del_adListShared",
	#"del_adListChunked",
	#"del_degAwareRHH",
	
	"threads2_graphite",
	"threads2_stinger",
	"threads2_adListShared",
	"threads2_adListChunked",
	
	"threads4_graphite",
	"threads4_stinger",
	"threads4_adListShared",
	"threads4_adListChunked",
	"threads4_degAwareRHH",
	
	"threads8_graphite",
	"threads8_stinger",
	"threads8_adListShared",
	"threads8_adListChunked",
	"threads8_degAwareRHH",
	]

def build_insert_csv(rootDir):
	with open(rootDir + "/update.csv", "w") as outf:
		for dataset in activeDatasetList:
			headerPrinted = False
			header = datasetMap[dataset] + ",batch,"
			for ds in activeDSList:
				fname = rootDir + '/bfsdyn/' + ds + '/' + dataset + '/Run3/Update.csv'
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
				for ds in activeDSList:
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Run3/Alg.csv'
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
			
			for ds in activeDSList:
				outf.write(algoMap[algo] + ',' + dsMap[ds])
			
				for dataset in activeDatasetList:
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Run3/Alg.csv'
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
			
			for ds in activeDSList:
				outf.write(algoMap[algo] + ',' + dsMap[ds])
			
				for dataset in activeDatasetList:
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Run3/Update.csv'
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
			
			for ds in activeDSList:
				outf.write(algoMap[algo] + ',' + dsMap[ds])
			
				for dataset in activeDatasetList:
					tot = 0.0
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Run3/Update.csv'
					print("Processing " + fname)
					with open(fname, "r") as csvf:				
						for line in csvf:
							tot = tot + float(line)
					fname = rootDir + '/' + algo + '/' + ds + '/' + dataset + '/Run3/Alg.csv'
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
#	build_summed_algo_csv(rootDir)
	build_summed_update_csv(rootDir)
#	build_summed_tot_csv(rootDir)	

	
	
		
		
