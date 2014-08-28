/*
 * Test.cpp
 *
 *  Created on: May 5, 2014
 *      Author: yqian33
 */

#include <iostream>
#include <fstream>
#include <string>
#include "def.h"
#include "auxilliary.h"
#include "networkblast.h"
#include "SemSim.h"
#include "anemo.h"
#include "netaligner.h"

using namespace std;

void testNetworkblast()
{
	cout<<"**************************************TEST NetworkBLAST*********************************\n";

	string dirprefix="/home/yqian33/APPI/verify_alignment_after_interaction_added/NetworkBLAST/yeast-worm-i2d/";
	string dir = dirprefix+"complexes/";
	string s2s="yeast-worm";
	string oufile= dirprefix+s2s+"sol.txt";

	if (mergeSolutionsNB(dir, oufile))
	{
		cout << "done" << endl;
	}

	map <int, string> complexMap = map <int, string> ();
	string infile = "/home/yqian33/APPI/software/AlignMCL-1.2/protein_complexes/yeast.txt";
	if (buildComplexMap(infile, complexMap))
	{
		cout << "Known complexes map constructed！" <<endl;
	}

	map <int, string> proMap = map <int, string> ();
	if (buildSolMap(oufile, proMap, 0))
	{
		cout << "Solutions map constructed！" <<endl;
	}

	string outfile0 = dirprefix + s2s + "_yeast_overlap_sol.txt";
	string outfile1 = dirprefix + s2s + "_yeast_overlap_complex.txt";
	string prefix=" ";
	if (overlap_Solution(proMap, complexMap, outfile0, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}
	if (overlap_Complex(proMap, complexMap, outfile1, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}

	//string solFile="/home/yqian33/APPI/software/networkblast/result/yeast-worm-MCL/yeast-worm-db.txt";
	/*string gofile="/home/yqian33/APPI/software/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_yeast";
	string outfile= dirprefix+"yeast0-worm_goterms.txt";

	if (compute_GOtermsPerProtein(oufile, gofile, outfile, 0))
	{
		cout << "done" << endl;
	}*/

	  /* string bpFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/go_terms.mgi.txt";
	string gotermsFile="/home/yqian33/APPI/software/networkblast/result/yeast-fly-MCL/1/yeast1-fly-dip-psi_goterms.txt";
	string out="/home/yqian33/APPI/software/networkblast/result/yeast-fly-MCL/1/yeast1-fly-dip-psi_goterms_bp.txt";

	if( pickBPtermsOfSolutions(gotermsFile, bpFile, out))
	{
		cout<<"done"<<endl;
	}*/
	cout<<"**************************************TEST NetworkBLAST*********************************"<<endl;
}

//COMPUTE intra specises semantic similarity for networkblast solutions
void testNetworkblastSS()
{
	string dirprefix="/home/yqian33/APPI/software/networkblast/result/yeast-worm-i2d-blastpdb/";
	string solFile= dirprefix +"yeast-worm-sol.txt";
	string dirname =dirprefix+"solutions/";
	if(computeIntraSpeciesSimGIC(solFile, 0, dirname))//computer intra species SS
	{
		cout<<"done"<<endl;
	}
}
//test semantic similarity for Biological Process
void testSSforBP()
{
	string solFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/AlignMCL/fly.DIP-yeast.DIP.txt";
	string dirname ="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/solutions/fly-yeast-DIP_fly/";
	if(computeIntraSpeciesSimGIC(solFile, 0, dirname))//computer intra species SS
	{
		cout<<"done"<<endl;
	}
}
/*
 * Annotate GO terms to proteins in solution of AlignMCL
 */
void testAlignMCL_GA()
{
	string solFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/AlignMCL/fly.DIP-yeast.DIP.txt";
	string gofile="/home/yqian33/APPI/software/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_fly";
	string outfile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/AlignMCL/fly-yeast-dip-fly_goterms.txt";

	if (compute_GOtermsPerProtein(solFile, gofile, outfile, 0))
	{
		cout << "done" << endl;
	}
}
void testAlignMCL()
{
	cout<<"**************************************TEST AlignMCL*********************************\n";
	map <int, string> complexMap = map <int, string> ();
	string infile = "/home/yqian33/APPI/software/AlignMCL-1.2/protein_complexes/fly.txt";

	map <int, string> complexMap2 = map <int, string> ();
	string infile2 = "/home/yqian33/APPI/software/AlignMCL-1.2/protein_complexes/yeast.txt";
	if (buildComplexMap(infile, complexMap))
	{
		cout << "Known complexes map constructed！" <<endl;
	}
	if (buildComplexMap(infile2, complexMap2))
	{
		cout << "Known complexes map constructed！" <<endl;
	}

	string dir="/home/yqian33/APPI/verify_alignment_after_interaction_added/AlignMCL/testremove/fly-yeast-rm10000/11/";
	map <int, string> proMap = map <int, string> ();
	string infile1=dir+"fly-yeast.txt";
	if (buildSolMap(infile1, proMap, 0))
	{
		cout << "Solutions map constructed！" <<endl;
	}

	map <int, string> proMap2 = map <int, string> ();
	//string infile1=dir+"fly-yeast.txt";
	if (buildSolMap(infile1, proMap2, 1))
	{
		cout << "Solutions map constructed！" <<endl;
	}

	string outfile =dir+"fly-yeast_fly_overlap_sol.txt";
	string outfile1=dir+"fly-yeast_fly_overlap_complex.txt";
	string prefix=" ";
	if (overlap_Solution(proMap, complexMap, outfile, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}
	if (overlap_Complex(proMap, complexMap, outfile1, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}

	string outfile_1 =dir+"fly-yeast_yeast_overlap_sol.txt";
	string outfile1_1=dir+"fly-yeast_yeast_overlap_complex.txt";
	//string prefix=" ";
	if (overlap_Solution(proMap2, complexMap2, outfile_1, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}
	if (overlap_Complex(proMap2, complexMap2, outfile1_1, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}



	cout<<"**************************************Done*********************************\n";
}
void testMAWISH()
{
	cout<<"**************************************TEST MaWish*********************************\n";

	string dir="/home/yqian33/APPI/verify_alignment_after_interaction_added/MaWISH/fly-yeast-DIP/";
	map <int, string> complexMap = map <int, string> ();
	string infile = "/home/yqian33/APPI/software/AlignMCL-1.2/protein_complexes/yeast.txt";
	if (buildComplexMap(infile, complexMap))
	{
		cout << "Known complexes map constructed！" <<endl;
	}

	map <int, string> proMap = map <int, string> ();
	string infile1=dir+"72271_49321.txt";
	if (buildSolMap(infile1, proMap, 1))
	{
		cout << "Solutions map constructed！" <<endl;
	}

	string outfile =dir+"fly-yeast_yeast_overlap_sol.txt";
	string outfile1=dir+"fly-yeast_yeast_overlap_comp.txt";
	string prefix="fly-yeast-i2d_yeast";
	if (overlap_Solution(proMap, complexMap, outfile, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}
	if (overlap_Complex(proMap, complexMap, outfile1, prefix, 3))
	{
		cout << "Finished!"<<endl;
	}

	//string gofile="/home/yqian33/APPI/ProSeqBla/idmapfile/DROME_7227_idmapping_selected.tab";
	/*string gofile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_yeast";
	string outfile3="/home/yqian33/APPI/software/MAWISH/data/alignment/"+dir+"/fly-yeast_yeast-goterms.txt";

	if (compute_GOtermsPerProtein(infile1, gofile, outfile3, 1))
	{
			cout << "done" << endl;
	}

	string dirname ="/home/yqian33/APPI/software/MAWISH/data/alignment/"+dir+"/solutions/";
	if(solution2queryList(infile1, 1, dirname))//computer intra species SS
	{
		cout<<"done"<<endl;
	}*/
	cout<<"**************************************TEST MaWish*********************************\n";

}

void testNetAligner()
{
	cout<<"**************************************TEST NetAligner*********************************\n";
	string s2s="yeast-fly";
	string dir="/home/yqian33/APPI/software/NetAligner_unix/test/"+s2s+"/result/";
	string oufile="/home/yqian33/APPI/software/NetAligner_unix/test/"+s2s+"/"+s2s+".sol";
	mergeSolutions(dir, oufile);

	map <int, string> complexMap = map <int, string> ();
	string infile = "/home/yqian33/APPI/software/AlignMCL-1.2/protein_complexes/yeast.txt";
	if (buildComplexMap(infile, complexMap))
	{
		cout << "Known complexes map constructed！" <<endl;
	}

	map <int, string> proMap = map <int, string> ();
	//string infile1="/home/yqian33/APPI/software/MAWISH/data/alignment/4932_7227_formatted.txt";
	if (buildSolMap(oufile, proMap, 0))
	{
		cout << "Solutions map constructed！" <<endl;
	}

	string outfile ="/home/yqian33/APPI/software/NetAligner_unix/test/"+s2s+"/Solution_overlap.txt";
	string outfile1="/home/yqian33/APPI/software/NetAligner_unix/test/"+s2s+"/Complex_overlap.txt";
	string prefix=s2s+"-i2d_yeast";
	if (overlap_Solution(proMap, complexMap, outfile, prefix, 4))
	{
		cout << "Finished!"<<endl;
	}
	if (overlap_Complex(proMap, complexMap, outfile1, prefix, 4))
	{
		cout << "Finished!"<<endl;
	}
	cout<<"**************************************TEST NetAligner*********************************\n";

}
int main(void)
{
	//testNetworkblast();
	//testMAWISH();
	testAlignMCL();


	//testNetworkblastSS();
	//testSSforBP();
	//testAlignMCL_GA();
	//test();
	//statistic_known( );
	//testNetAligner();

	return 0;
}


/*-------------------------------------------------------------------
 * --------------------------Merge solutions of Networkblast---------*/

	/*string dir = string("/home/yqian33/APPI/software/networkblast/result/fly-yeast-Nemo/complexes/");
	string oufile="/home/yqian33/APPI/software/networkblast/result/fly-yeast-Nemo/yeast-fly.txt";
	if (mergeSolutionsNB(dir, oufile))
	{
		cout << "done" << endl;
	}*/



/*-------------------------------------------------------------------------
 * ------------------------Find overlap with known complexes---------------
 */

	/*map <int, string> complexMap = map <int, string> ();
	string infile = "/home/yqian33/APPI/software/AlignMCL-1.2/protein_complexes/yeast.txt";
	if (buildComplexMap(infile, complexMap))
	{
		cout << "Known complexes map constructed！" <<endl;
	}

	map <int, string> proMap = map <int, string> ();
	string infile1="/home/yqian33/APPI/software/MAWISH/data/alignment/4932_7227_formatted.txt";
	if (buildSolMap(infile1, proMap, 0))
	{
		cout << "Solutions map constructed！" <<endl;
	}

	string outfile ="/home/yqian33/APPI/software/MAWISH/data/alignment/yeast-fly_yeast_overlap_sol.txt";
	string outfile1="/home/yqian33/APPI/software/MAWISH/data/alignment/yeast-fly_yeast_overlap_complex.txt";
	string prefix="yeast.DIP-fly.DIP_yeast";
	if (overlap_Solution(proMap, complexMap, outfile, prefix, 4))
	{
		cout << "Finished!"<<endl;
	}
	if (overlap_Complex(proMap, complexMap, outfile1, prefix, 4))
	{
		cout << "Finished!"<<endl;
	}*/

	//statistic_known();



/*-------------------------------------------------------------------
 * -------------------for networkblast test--------------------------
 */
	//string infile="/home/yqian33/APPI/software/networkblast/MIPSyeast_complex.txt";
	//string outfile="/home/yqian33/APPI/software/networkblast/MIPSyeast_complex_clean.txt";
	//map <string, string> complexMap = map <string, string> ();
	//mipsComplexProcess(infile, outfile, complexMap);

	/*string overlapFile="/home/yqian33/APPI/software/networkblast/result/yeast-fly-Nemo/yeast-fly.txt";
	string outfile="/home/yqian33/APPI/software/networkblast/result/yeast-fly-Nemo/pureComplex-yeast-fly_yeast.txt";
	string complexFile="/home/yqian33/APPI/software/networkblast/MIPSyeast_complex.txt";
	if( pureComplexVerify(overlapFile, complexFile, outfile, 0) )
	{
		 cout <<"done" << endl;
	}*/


/*-------------------------------------------------------------------
 * -----------------------GO terms annotation test-------------------
 * find all the GO terms annotated to every proteins
 */
	/*string solFile="/home/yqian33/APPI/software/networkblast/result/fly-yeast-Nemo/fly-yeast.txt";
	string gofile="/home/yqian33/APPI/software/fastSemSim-0.8.1/test/data/ACs/gene_association.goa_yeast";
	string outfile="/home/yqian33/APPI/software/networkblast/result/fly-yeast-Nemo/fly-yeast-dip-yeast_goterms.txt";

	if (compute_GOtermsPerProtein(solFile, gofile, outfile, 1))
	{
		cout << "done" << endl;
	}*/

	/*string solFile="/home/yqian33/APPI/software/MAWISH/data/alignment/4932_6239.txt";
	string gofile="/home/yqian33/APPI/ProSeqBla/idmapfile/YEAST_559292_idmapping_selected.tab";
	string outfile="/home/yqian33/APPI/software/MAWISH/data/alignment/4932_6239-yeast-goterms.txt";

	if (compute_Goterms_MaWISh(solFile, gofile, outfile, 0))
	{
			cout << "done" << endl;
	}*/



/*--------------------------------------------------------------------------------------
 * -----------------------Semantic Similarity test  intra and interSS-------------------
 */

	/*string solFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/AlignMCL/yeast.DIP-worm.DIP.txt";
	string dirname ="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/solutions/yeast.DIP-worm.DIP-yeast/";
	if(solution2queryList(solFile, 0, dirname))//computer intra species SS
	{
		cout<<"done"<<endl;
	}*/

	/*string solFile="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/AlignMCL/fly.DIP-yeast.DIP.txt";
	string dirname ="/home/yqian33/APPI/software/AlignMCL-1.2/alignments/solutions/fly.DIP-yeast.DIP/";
	if( computeInterSpeciesSimGIC(solFile, dirname) )
	{
		cout<<"done"<<endl;
	}*/

	//testAnemo();
	//test();


